from typing import Iterable, Sequence, Tuple

import numpy as np
from cirq.ops import gate_features
import cirq

NOISY_OP_TAG = "e"

@cirq.value.value_equality
class AsymmetricBitflip(gate_features.SingleQubitGate):
    """Flip bits with unequal probabilities.

    """

    def __init__(self, p0: float, p1: float) -> None:
        r"""The asymmetric classical bit flip channel.

        Construct a channel to model simultaneous relaxation and excitation
        according to different probabilities. This is different from a
        composition of relaxation/excitation channels, and can furthermore
        not be modelled by a composition of symmetric bitflip and amplitude
        damping channels.

        This channel evolves a density matrix via

            $$
            \rho \rightarrow M_0 \rho M_0^\dagger
                           + M_1 \rho M_1^\dagger
                           + M_2 \rho M_2^\dagger
            $$

        With

            M_0 = \sqrt{p0} |1><0|
            M_1 = \sqrt{p1} |0><1|
            M_2 =& \sqrt{1-p0} * |0><0| + \sqrt{1-p1}|1><1|


        Args:
            p0: the probability of |0> -> |1>
            p1: the probability of |1> -> |0>

        Raises:
            ValueError: if gamma or p is not a valid probability.
        """
        self._p0 = cirq.value.validate_probability(p0, 'p0')
        self._p1 = cirq.value.validate_probability(p1, 'p1')

    def _channel_(self) -> Iterable[np.ndarray]:
        sqrt_p0 = np.sqrt(self._p0)
        sqrt_p1 = np.sqrt(self._p1)
        M0 = np.array([[0., 0.], [sqrt_p0, 0.]])
        M1 = np.array([[0., sqrt_p1], [0., 0.]])
        M2 = np.asarray([[np.sqrt(1-self._p0), 0], [0, np.sqrt(1-self._p1)]])
        assert np.allclose(M0.T @ M0 + M1.T @ M1 + M2.T @ M2, np.eye(2))
        return (M0, M1, M2)

    def _has_channel_(self) -> bool:
        return True

    def _value_equality_values_(self):
        return self._p0, self._p1

    def __repr__(self) -> str:
        return 'cirq.asymmetric_bitflip(p0={:1.3f},p1={:1.3f})'.format(
            self._p0, self._p1
        )

    def __str__(self) -> str:
        return 'asymmetric_bitflip(p0={!r},p1={!r})'.format(self._p0, self._p1)

    def _circuit_diagram_info_(self, args: cirq.protocols.CircuitDiagramInfoArgs
    ) -> str:
        if args.precision is not None:
            f = '{:.' + str(args.precision) + 'g}'
            return 'ABF({},{})'.format(f, f).format(self._p0, self._p1)
        return 'ABF({!r},{!r})'.format(self._p0, self._p1)


def asymmetric_bitflip(p0: float, p1: float) -> AsymmetricBitflip:
    return AsymmetricBitflip(p0=p0, p1=p1)


@cirq.value.value_equality
class TwoQubitSymmetricDepolarizingChannel(cirq.ops.gate_features.TwoQubitGate):
    """
    A channel that depolarizes symmetrically for two-qubit.
    """

    def __init__(self, p: float) -> None:
        """
        Construct a depolarizing channel with p/16 probability of depolarization
        along any 2-qubit Pauli term.

        Args:
            p: Probability of incoherent noise.

        """
        self._p = cirq.value.validate_probability(p, 'p')
        # probability of any given 2-term Pauli noise
        self._p_op = p / 16.
        # probability of Identity
        self._p_i = 1 - 15. * p / 16.

    def _mixture_(self) -> Sequence[Tuple[float, np.ndarray]]:
        """Define the mixture of all two-qubit pauli terms.

        This mixture defines 16 Kraus operators M_i, taken from the set

            M_i = {(1-p*15/16)*I, 1/16*IX, 1/16 IY, 1/16 IZ, 1/16 XX, ...}

        And the result of applying the channel sends rho to

            M_0 rho M_0.T + M_1 rho M_1.T + ...
        """
        q = cirq.LineQubit(0) # dummy qubit
        pauli_gates = [cirq.X, cirq.Y, cirq.Z]
        pauli_strings = [cirq.PauliString({q: p}) for p in pauli_gates]
        paulis = [np.eye(2, dtype=np.complex256)] + [cirq.unitary(p).astype(np.complex256) for p in pauli_strings]
        # paulis = [np.eye(2)] + [cirq.unitary(p) for p in pauli_strings]
        double_paulis = [np.kron(m1, m2).reshape(4, 4) for m1 in paulis for m2 in paulis]
        mixture_tuple = [(self._p_i, double_paulis[0])]
        for op in double_paulis[1:]:
            mixture_tuple += [(self._p_op, op)]
        mixture_tuple = tuple(mixture_tuple)
        atol = 1e-2
        rtol = 1e-2
        lhs = np.sum(m[1].T.conj() @ m[1] * m[0] for m in mixture_tuple)
        rhs = np.eye(4)
        if not np.allclose(lhs, rhs, rtol=rtol, atol=atol):
            diff = np.sum(rhs - lhs)
            print(f" WARNING: channel not within {rtol} rtol, {atol} atol of trace preserving. \
                    Net matrix difference = {diff}")
            # raise ValueError(f"channel not within {rtol} rtol, {atol} atol of trace preserving. \
            #         Net matrix difference = {diff}")
        return mixture_tuple

    def _has_mixture_(self) -> bool:
        return True

    def _value_equality_values_(self):
        return self._p

    def __repr__(self) -> str:
        return 'two_qubit_symmetric_depolarize(p={:1.3f})'.format(self._p)

    def __str__(self) -> str:
        return self.__repr__()

    def _circuit_diagram_info_(
            self, args: cirq.protocols.CircuitDiagramInfoArgs) -> str:
        if args.precision is not None:
            f = '{:.' + str(args.precision) + 'g}'
            gate_msg = '2Q-S-D({})'.format(f).format(self._p)
        else:
            gate_msg = '2Q-S-D({})'.format(self._p)

        output = cirq.protocols.CircuitDiagramInfo(
            wire_symbols=(gate_msg, gate_msg))
        return output


def symmetric_depolarize(p: float) -> TwoQubitSymmetricDepolarizingChannel:
    return TwoQubitSymmetricDepolarizingChannel(p)


class TwoQubitGateNoiseModel(cirq.devices.NoiseModel):
    """Applies noise to each qubit individually at the end of every moment."""

    def __init__(self, qubit_noise_channel: 'cirq.Channel'):
        # if qubit_noise_channel.num_qubits() != 2:
        #     raise ValueError('noise.num_qubits() != 2')
        self.qubit_noise_channel = qubit_noise_channel

    def noisy_operation(self, *operation: 'cirq.Operation'):
        return self.qubit_noise_channel(*operation)


def DephasingNoiseWrapper(circuit: cirq.Circuit, gamma: float) -> cirq.Circuit:
    """Wrap an existing circuit with the `two_qubit_phase_damp` noise model.

    This will copy an existing circuit and apply the custom two-qubit gate
    dephasing noise specified in TwoQubitPhaseDampingChannel to all ZZPowGate
    objects in an existing circuit, returning the new circuit.

    Args:
        circuit (cirq.Circuit): Existing circuit to wrap with noise (this is
            not modified!)
        gamma: The damping constant.
    """
    noise_channel = two_qubit_phase_damp(gamma)
    noise_model = TwoQubitGateNoiseModel(noise_channel)

    noisy_circuit = cirq.Circuit()
    for moment in circuit:
        # carry over original circuit gates
        noisy_circuit.append(moment)
        # tack on noise to ZZ gates
        for op in moment.operations:
            # FIXME: Can't access gate identity or superclass; op.__bases__
            # isn't seeing the inheritance on ZZPowGate...
            if len(op.qubits) == 2:
                noisy_circuit.append(noise_model.noisy_operation(*op.qubits))

    return noisy_circuit


def apply_noise(circuit: cirq.Circuit,
                noise_model: cirq.NoiseModel) -> cirq.Circuit:
    """A placeholder for applying a noise model in Cirq.

    For a pairwise kernel, it is assumed that there is some optimization scheme
    for "dense-packing"  gates in the kernel ansatz (i.e. every qubit interacts
    with a gate at  every moment) and therefore single-qubit noise is only
    applied after each gate instead of at every moment. This results in the same
    total amount of dephasing as if the gates had been optimally packed.
    """

    noisy_circuit = cirq.Circuit()
    for moment in circuit:
        # carry over original circuit gates, plus noisy moment
        noisy_circuit.append(noise_model.noisy_moment(moment, moment.qubits))

    return noisy_circuit


def random_2q_pauli_exp(p: float, qubits: Iterable[cirq.Qid]) -> cirq.Operation:
    """Generate a random 2-qubit Pauli of the form exp(i*p*P).

    P is sampled from the set of all 2-qubit paulis, minus the Identity. i.e.

        {P} = {X0, Y0, Z0, X0*X1, Y0*X1, ... Z0*Z1} \ {I}

    where '*' denotes tensor product.
    """
    P = cirq.PauliString(cirq.I(qubits[0]))
    while np.allclose(cirq.unitary(P), 1):
        P0 = np.random.choice([cirq.I, cirq.X, cirq.Y, cirq.Z]).on(qubits[0])
        P1 = np.random.choice([cirq.I, cirq.X, cirq.Y, cirq.Z]).on(qubits[1])
        P = cirq.PauliString(P0) * cirq.PauliString(P1)
    return np.exp(1j * p * P)


def decay_prob(T1, T2, gate_time):
    """Determine the decay prob from T1 and T2 in gate_time.

    Given a gate execution duration and device T1, T2 times, compute a
    probability for decay under an amplitude damping channel (corresponding to
    T1) and under a dephasing channel (corresponding to T2).

    WARNING: This will clash with any channel simulating noise corresponding to
    a single-qubit gate fidelity, as any diagnostic that determins gate fidelity
    has implicitly accounted for time decay. Loosely,

        (observed gate infidelity) = (execution time decay) + (gate infidelity)

    To get around this clash, a noise model should implement logic such that
    any channels that are parameterized based on single qubit gate fidelity
    will overwrite channels acting on the same operation based on T1,T2 decay.
    """
    p_decay = 0
    if not np.isclose(T1, 0):
        p_decay += gate_time / T1

    # The T2 quoted from experiment actually includes T1 AND T2 decay; therefore
    # the decay probability needs to be reduced and used to solve for just
    # decay prob due to T2
    p_dephasing = 0
    if not np.isclose(T2, 0):
        p_dephasing = (gate_time / T2) / (1 - 4/3 * p_decay)

    return(p_decay, p_dephasing)