import itertools
from typing import Sequence, Optional

import numpy as np
import cirq

import simple_noise_model
import device_specs
import channels


class TwoLocalDepolarizeAndDephase(cirq.devices.NoiseModel):
    """A two-local depolarizing noise model with ancilla dephasing/amplitude damp

    Given parameters:
        T1: Ancilla T1 (us) relaxation time.
        T2_eff: Effective T2 (us) for ancilla, after any dynamical decoupling scheme
        err_1q: Characteristic error probability for single-qubit gates
        err_2q: Characteristic error probability for two-qubit gates

    Perform the following
        1. For every single qubit gate on every qubit: append a SINGLE
            QUBIT DEPOLARIZING CHANNEL parameterized by prob ~ `err_1q`.
        2. For every two qubit operation, apply a SYMMETRIC TWO QUBIT
            DEPOLARIZING CHANNEL parameterized by parameterized by prob ~ `err_2q`.
        3. For every circuit timestep dt:
            - apply a dephasing channel on the ancilla, parameterized by
                `p_dephasing * dt`
            - apply a amplitude damping channel on the ancilla, parameterized by
                `p_decay * dt`
            where `(p_dephasing, p_decay)` are returned according to the logic
                of `channels.decay_prob`
        """
    def __init__(self, ancilla: cirq.Qid, err_1q: float, err_2q: float, T1: float, T2_eff: float, sampled_error: Optional[bool] = True, verbose: Optional[bool]=True):
        """Generate a generic error model that depolarizes one- and two-qubit gates symmetrically.

        Args:
            ancilla: Indicates which qubit is the Ancilla
            err_1q: characteristic error for 1q depol channels
            err_2q: characteristic error for 2q depol channels
            T1, T2_eff: T1 and effective T2 (after dynamic decoupling) for
                ancilla. Note that a value of 0 may be used to indicate "no error".
            sampled_error: Sample each parameter as err*Uniform[0,1]. If false, every single noise channel
                will have parameter exactly err_1q or err_2q for one- and two-qubit gates respectively.
        """
        self.ancilla = ancilla
        self.err_1q = err_1q
        self.err_2q = err_2q
        self.T1 = T1
        self.T2_eff = T2_eff
        self.sampled_error = sampled_error
        self.verbose = verbose

    def _make_twolocal_noise_model(self, circuit):
        # Prepare a template noise model from the single/two qubit errors provided
        qubits = list(circuit.all_qubits())
        nqubits = len(qubits)

        qubit_keys = [tuple([q]) for q in qubits]
        if self.sampled_error:
            err_1q_lst = np.random.random(nqubits) * self.err_1q
        else:
            err_1q_lst = np.ones(nqubits) * self.err_1q

        err_1q_map = dict(zip(qubit_keys, err_1q_lst))
        # Some two-qubit gate errors
        qubit_pairs = list(itertools.combinations(qubits, 2))# all-to-all connectivity
        # FIXME: I'm being lazy here and not enforcing consistent edge ordering.
        qubit_pairs += [(edge[1], edge[0]) for edge in qubit_pairs]

        if self.sampled_error:
            err_2q_lst = np.random.random(len(qubit_pairs)) * self.err_2q
        else:
            err_2q_lst = np.ones(len(qubit_pairs)) * self.err_2q
        err_2q_map = dict(zip(qubit_pairs, err_2q_lst))

        # I'll completely ignore two-qubit unitary error in this model
        calibration_dct = {"fid_1q": err_1q_map, "fid_2q_pure": err_2q_map, "fid_2q_tot": err_2q_map}
        return simple_noise_model.SimpleNoiseModel.from_calibration_data(calibration_dct, verbose=self.verbose)

    def noisy_moments(self, circuit: cirq.Circuit,
                      system_qubits: Sequence[cirq.Qid]
                     ) -> Sequence[cirq.OP_TREE]:
        """Construct a two-local depol model with ancilla damping/dephasing.

        We make a very rough approximation that the circuit duration (ns) is:
            4000 + 25 * depth
        """
        twolocal_model = self._make_twolocal_noise_model(circuit)
        noisy_circuit_first_pass = circuit.with_noise(twolocal_model)

        # Instead of many timesteps of channels, just make as single channel
        # with decay probability given by the entire circuit
        circuit_duration = 0
        for moment in circuit:
            circuit_duration += max(device_specs.get_gate_duration(gate, verbose=self.verbose) for gate in moment)
        # Convert to us
        circuit_duration *= 1e-3
        # circuit_duration += 4000 # measurement duration, ns

        new_op_tree = []
        for i, moment in enumerate(noisy_circuit_first_pass):
            new_op_tree.append(moment)
            # Sort of assume any gate preparation circuit plus stacked noise ops
            # will not make it this deep
            if i == 5:
                p_decay, p_dephasing = channels.decay_prob(self.T1, self.T2_eff, circuit_duration)
                if not np.isclose(self.T1, 0):
                    new_op_tree.append(cirq.amplitude_damp(gamma=p_decay).on(self.ancilla))
                if not np.isclose(self.T2_eff, 0):
                    new_op_tree.append(cirq.phase_damp(gamma=p_dephasing).on(self.ancilla))

        return new_op_tree