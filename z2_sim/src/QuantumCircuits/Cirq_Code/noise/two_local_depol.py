import itertools
from typing import Sequence, Optional

import numpy as np
import cirq

from z2_sim.src.QuantumCircuits.Cirq_Code.noise import simple_noise_model
from z2_sim.src.QuantumCircuits.Cirq_Code.noise import channels

class TwoLocalDepol(cirq.devices.NoiseModel):
    """A two-local depolarizing noise model.

    Given parameters:
        err_1q: Characteristic error probability for single-qubit gates
        err_2q: Characteristic error probability for two-qubit gates

    Perform the following
        1. For every single qubit gate on every qubit: append a SINGLE
            QUBIT DEPOLARIZING CHANNEL parameterized by prob ~ `err_1q`.
        2. For every two qubit operation, apply a SYMMETRIC TWO QUBIT
            DEPOLARIZING CHANNEL parameterized by parameterized by prob ~ `err_2q`.
        """
    def __init__(self, err_1q: float, err_2q: float, sampled_error: Optional[bool] = True, verbose: Optional[bool]=True):
        """Generate a generic error model that depolarizes one- and two-qubit gates symmetrically.

        Args:
            err_1q: characteristic error for 1q depol channels
            err_2q: characteristic error for 2q depol channels
        """
        self.err_1q = err_1q
        self.err_2q = err_2q
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
        """Construct a two-local depol model."""


        out = []
        for moment in circuit:
            new_moment = []

            for gate in moment:
                if channels.NOISY_OP_TAG in gate.tags:
                    continue
                gate_qubits = tuple(gate.qubits)
                n_qubit = len(gate_qubits)

                if n_qubit == 1:
                    tagged_op = cirq.TaggedOperation(cirq.depolarize(p=self.err_1q).on(gate_qubits[0]), channels.NOISY_OP_TAG)
                elif n_qubit == 2:
                    tagged_op = cirq.TaggedOperation(channels.symmetric_depolarize(p=self.err_2q).on(*gate_qubits), channels.NOISY_OP_TAG)

                new_moment.append(tagged_op)
            new_op_tree = [moment, new_moment]
            out.append(new_op_tree)

        return out

        # twolocal_model = self._make_twolocal_noise_model(circuit)
        # noisy_circuit_first_pass = circuit.with_noise(twolocal_model)
        #
        #
        #
        # # OP_TREE signature, not sure if necessary tbh
        # return [moment for moment in noisy_circuit_first_pass]
