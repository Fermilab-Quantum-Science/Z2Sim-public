import numbers
import itertools

from typing import Sequence, Optional

import numpy as np
import cirq


from z2_sim.src.QuantumCircuits.Cirq_Code.noise import channels


class ZZCrossTalkV1(cirq.devices.NoiseModel):
    """A basic crosstalk simulation

    Given parameters:
        zeta_dict: characteristic strength for cross-Kerr interactions between
            qubits (i, j)

    Perform the following
        1. For every two qubit operation, apply the ZZ interaction gate

            exp( (-i2π ζ_{ij} T )* |11><11|  )

         parameterized by parameterized by `Zeta_{ij}` and gate duration T.
        """
    TARGET_GATES_IMPLEMENTED = [
    "CNOT",
    "SIS"
    ]
    def __init__(self, zeta, target_gate, gate_duration=0.01, sampled_error: Optional[bool] = False, verbose: Optional[bool]=True):
        """Generate a generic error model that depolarizes one- and two-qubit gates symmetrically.

        Args:
            zeta: Characterization of crosstalk strength. This can either be a
                constant (in which case all pairs of qubits will be subjected
                to the same noise model) or a dictionary maping qubit pairs to
                floats.
            target_gate: A string describing which two-qubit gate to
                apply the crosstalk noise to. Currently supported:
                    `CNOT`
            gate_duration: Typical gate execution time in microseconds.
            sampled_error: Sample each parameter as err*Uniform[0,1]. If false, every single noise channel
                will have parameter exactly err_1q or err_2q for one- and two-qubit gates respectively.
        """

        self.zeta = zeta
        self.zeta_map = None

        if target_gate not in self.TARGET_GATES_IMPLEMENTED:
            raise NotImplementedError(
                "Noise model currently does not support "
                "target gate: {}".format(target_gate)
            )
        self.target_gate = target_gate
        self.gate_duration = gate_duration
        self.sampled_error = sampled_error
        self.verbose = verbose

    def _make_zeta_map(self, zeta: float, circuit: cirq.Circuit):
        """Construct a map from a flat noise parameter."""

        if isinstance(zeta, dict):
            return zeta

        if isinstance(zeta, numbers.Number):
            qubits = list(circuit.all_qubits())

            qubit_pairs = list(itertools.combinations(qubits, 2))# all-to-all connectivity
            qubit_pairs += [(edge[1], edge[0]) for edge in qubit_pairs]

            if self.sampled_error:
                err_2q_lst = np.random.random(len(qubit_pairs)) * zeta
            else:
                err_2q_lst = np.ones(len(qubit_pairs)) * zeta
            return dict(zip(qubit_pairs, err_2q_lst))

        raise ValueError("Noise model must be initialized with a `zeta` that is"
                         " either numeric or a map of qubit pairs to numerics.")

    def _is_target_gate(self, gate: cirq.Gate):
        """Construct a predicate for `gate` based on the target gate."""

        if self.target_gate == "CNOT":
            return (isinstance(gate, cirq.CXPowGate) and gate.exponent == 1)
        elif self.target_gate == "SIS":
            return (isinstance(gate, cirq.ISwapPowGate) and gate.exponent == 0.5)

        return False

    def noisy_moments(self, circuit: cirq.Circuit,
                      system_qubits: Sequence[cirq.Qid]
                     ) -> Sequence[cirq.OP_TREE]:
        """Instantiate the noise model by modifying an existing circuit.

        Args:
            circuit: The circuit to add noise to.
            system_qubits: For callback compatibility with `Circuit.with_noise`

        Returns:
            A sequence of OP_TREEs enacting all noise channels specified by
                nonzero input noise parameters on the input circuit.
        """
        # Delay this processing step until call time, so that we can access the
        # qubits in the circuit.
        self.zeta_map = self._make_zeta_map(self.zeta, circuit)
        out = []

        for i, moment in enumerate(circuit):
            # Attempting to get a clean separation between original gates
            # and new ops effecting noise
            new_moment = []

            # moment_duration = max(device_specs.get_gate_duration(gate, verbose=self.verbose) for gate in moment)
            for op in moment:
                gate_qubits = tuple(op.qubits)
                n_qubit = len(gate_qubits)
                gate = op.gate
                if (channels.NOISY_OP_TAG not in op.tags) and (n_qubit == 2) and self._is_target_gate(gate):
                    sorted_2q = tuple(sorted(gate_qubits, key=lambda q: q.row))
                    zeta = self.zeta_map.get(sorted_2q)
                    if zeta is None:
                        raise ValueError("Did not receive zeta value for qubits {}".format(sorted_2q))

                    if zeta > 0:
                        # We enforce tagging on noise operations to avoid
                        # undesirable "stacking" of independent noise models
                        tagged_op = cirq.TaggedOperation(cirq.cphase(rads=2*np.pi*self.gate_duration*zeta).on(*gate_qubits), channels.NOISY_OP_TAG)
                        new_moment.append(tagged_op)

            out.append([new_moment, moment])

        return out

