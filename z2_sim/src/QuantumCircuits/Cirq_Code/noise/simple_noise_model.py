from typing import Dict, Sequence, Tuple
import numpy as np
import cirq

from z2_sim.src.QuantumCircuits.Cirq_Code.noise import channels, device_specs


one_keys = [
    'p0_readout_error',
    'p1_readout_error',
    'single_qubit_readout_separation_error',
    'single_qubit_t1',
    'single_qubit_rb_total_error']
two_keys = [
    'fid_2q_tot',
    'fid_2q_pure']


class ConstantMap:
    """A very, very unsafe constant-valued dictionary."""
    def __init__(self, val):
        self.val = val

    def __getitem__(self, arg):
        return self.val

    def get(self, arg):
        return self.val



class SimpleNoiseModel(cirq.devices.NoiseModel):
    """Model a set of decoherence metrics via noise channels.

    Here is a procedural definition of this noise model:
        1. For every single qubit gate on every qubit:
            If a single qubit gate fidelity is less than 1, apply a SINGLE
            QUBIT DEPOLARIZING CHANNEL parameterized by the gate fidelity.
            Else, apply a SINGLE QUBIT DEPOLARIZING CHANNEL and a SINGLE QUBIT
            DEPHASING CHANNEL parameterized by the qubit T1, T2, and the longest
            hardware gate duration occuring in the same moment.
        2. For every empty timestep on every qubit:
            If T1 decay time is provided, apply a SINGLE QUBIT AMPLITUDE DAMPING
             CHANNEL parameterized by T1 and the longest hardware gate duration
             occuring in the same moment.
            If T2 decay time is provided, apply a SINGLE QUBIT DEPHASING CHANNEL
             parameterized by T1 and the longest hardware gate duration occuring
             in the same moment.
        3. For every two qubit operation:
            If a two qubit gate fidelity based on "pure" or control error is
            provided, apply a RANDOM TWO QUBIT ROTATION exp(iP_i * P_j)
            parameterized by the (tot_infidelity - pure_infidelity).
            If a two qubit gate fidelity containing total error is provided,
            provided, apply a SYMMETRIC TWO QUBIT DEPOLARIZING CHANNEL
            parameterized by the pure_infidelity.
            Else, apply a SINGLE QUBIT DEPOLARIZING CHANNEL and a SINGLE QUBIT
            DEPHASING CHANNEL parameterized by the qubit T1, T2, and the longest
            hardware gate duration occuring in the same moment.
        4. For every single qubit measurement:
            If measurement error rates are provided, apply an ASYMMETRIC BITFLIP
            CHANNEL parameterized by the readout error rates (p0, p1)
            Else, do nothing to the measurement.

    A set of high level recommendations:
        - p0/p1 errors are efficiently simulable within this framework but cannot
        account for general correlations in measurements. It might be preferable
        to disable p0/p1 errors and instead use "classical" postprocessing to
        apply a full (2 ** n, 2 ** n) response matrix to the output computational
        basis state probability distribution or observable.
        - The model for two-qubit gate error is incredibly sketchy. I recommend
        a simpler implementation, or none at all...
        - Single qubit gate fidelities are high enough that they might be
        well modeled by time decay alone.

    A set of significant limitations to the model
        1. No crosstalk exists in this model.
        2. I have not enabled gate-specific fidelities; this is highly unrealistic.
        3. Idle T1,T2 do not correspond to real T1, T2 - see (1)
        4. Symmetric depolarization isn't a realistic decoherence mechanism.

    A set of known bugs:
        1. FIXME: This is unlikely to be compatible with `cirq.measure_each`.
            It currently expects to find a distinct measurement operation for
            every qubit in the terminal moment.
    """

    def __init__(self, T1_map, T2_map,
                 p0_map, p1_map,
                 fid_1q_map,
                 fid_2q_tot_map, fid_2q_pure_map, verbose=True):
        """Initialize this noise model.

        All time-related quantities are in microseconds.

        Args:
            A set of maps that have Iterable[Qid] as keys, for example
                T1_map = {(cirq.LineQubit(0), ): 0.023}
                fid_2q_tot_map = {(cirq.LineQubit(0), cirq.LineQubit(1)): 0.97}
        """
        # If a dictionary is not provided in one of the constructors, it is
        # passed to init as `None` and gets gets converted to ConstantMap[0]
        self.T1_map = T1_map or ConstantMap(0.0)
        self.T2_map = T2_map or ConstantMap(0.0)
        # Probability p(0|1) for computational basis PVM {|0><0|, |1><1|}
        # for every qubit, individually
        self.p0_map = p0_map or ConstantMap(0.0)
        # Probability p(1|0)
        self.p1_map = p1_map or ConstantMap(0.0)
        # Single qubit gate infidelities
        self.fid_1q_map = fid_1q_map or ConstantMap(0.0)

        self.fid_2q_tot_map = fid_2q_tot_map or ConstantMap(0.0)
        self.fid_2q_pure_map = fid_2q_pure_map or ConstantMap(0.0)

        # Control the program whining about missing information
        self.verbose = verbose

    @classmethod
    def from_calibration_data(cls, calibration_dict: Dict[str, Dict[Tuple[cirq.Qid, ...], float]], verbose=True):
        """Initialize a noise model from a dict of calibration data.

        This dict can be generated from `load_calibration_data`. Currently
        this expects the following keys:

            'p0_readout_error' - Measurement probability of flipping '0' -> '1'
            'p1_readout_error' - Measurement probability of flipping '1' -> '0'
            'single_qubit_t1' - (Idle) qubit T1 (microseconds).
            'single_qubit_t2' - (Idle) qubit T2 (microseconds).
            'fid_1q' - Single qubit gate fidelities
            'fid_2q_tot' - two-qubit symmetric depolarizing channel parameter.
            'fid_2q_pure' - two qubit gate infidelity due to _incoherent_ noise.
                (I'm sorry for the legacy naming system, wasn't my idea.)

        """
        # TODO: validate input dict.
        T1_map = calibration_dict.get('single_qubit_t1')
        T2_map = calibration_dict.get('single_qubit_t2')
        if T2_map is None and verbose:
            print("WARNING: did not find any T2 data. Did you update "
                  "calibration data to include t2 values?")

        fid_2q_pure_map = calibration_dict.get(
            'fid_2q_pure')

        fid_1q_total_map = calibration_dict.get('fid_1q')
        fid_2q_tot_map = calibration_dict.get(
            'fid_2q_tot')
        p0_map = calibration_dict.get('p0_readout_error')
        p1_map = calibration_dict.get('p1_readout_error')

        return cls(T1_map, T2_map, p0_map, p1_map,
                   fid_1q_total_map,
                   fid_2q_tot_map=fid_2q_tot_map,
                   fid_2q_pure_map=fid_2q_pure_map, verbose=verbose)

    @classmethod
    def from_fixed_values(cls, T1=0.0, T2=0.0, p0=0.0, p1=0.0,
                          fid_1q=0.0, fid_2q_tot=0.0, fid_2q_pure=0.0, verbose=True):
        """Initialize a noise model from a set of constant parameters.

        If the quantity describes qubit noise properties, all qubits or qubit
        pairs in the circuit will be assigned that single quantity.

        Args:
            T1: (Idle) qubit T1 (microseconds).
            T2: (Idle) qubit T2 (microseconds).
            p0: Measurement probability of flipping '0' -> '1'.
            p1: Measurement probability of flipping '1' -> '0'.
            fid_1q: TODO
            fid_2q_tot: Total Sycamore Gate XEB infidelity.
            fid_2q_pure: Sycamore Gate XEB infidelity due to coherent noise.

        """
        T1_map = ConstantMap(T1)
        T2_map = ConstantMap(T2)
        p0_map = ConstantMap(p0)
        p1_map = ConstantMap(p1)
        fid_1q_map = ConstantMap(fid_1q)
        fid_2q_tot_map = ConstantMap(fid_2q_tot)
        fid_2q_pure_map = ConstantMap(fid_2q_pure)
        return cls(T1_map, T2_map, p0_map, p1_map, fid_1q_map,
                   fid_2q_tot_map=fid_2q_tot_map,
                   fid_2q_pure_map=fid_2q_pure_map, verbose=verbose)

    def __repr__(self):
        out = "noise model specs:"
        attrs_to_dump = [
            (self.T1_map, 'T1'),
            (self.T2_map, 'T2'),
            (self.p0_map, 'p0'),
            (self.p1_map, 'p1'),
            (self.fid_1q_map, 'fid_1q'),
            (self.fid_2q_tot_map, 'fid_2q_tot'),
            (self.fid_2q_pure_map, 'fid_2q_pure'),
        ]
        for attr, name in attrs_to_dump:
            if isinstance(attr, ConstantMap):
                out += "\n{}: {} (constant)".format(name, attr.get(None))
            elif isinstance(attr, dict):
                out += "\n{}: {}".format(name, attr)
            elif attr is None:
                out += "\n{}: None".format(name)
        return out

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
        out = []

        # TODO: validator to make sure the circuit contains only valid gates
        # for `get_gate_duration`
        last_moment = False
        for moment in circuit:
            # Attempting to get a clean separation between original gates
            # and new ops effecting noise
            new_moment = []
            # Any gates in this moment _don't_ get T1/T2 depolarized
            touched = set()

            moment_duration = max(device_specs.get_gate_duration(gate, verbose=self.verbose) for gate in moment)
            for gate in moment:
                if channels.NOISY_OP_TAG in gate.tags:
                    continue
                gate_qubits = tuple(gate.qubits)
                n_qubit = len(gate_qubits)
                # Assign measurement gates bit flip prob corresponding to
                # readout error probabilities.
                if isinstance(gate.gate, cirq.MeasurementGate):
                    last_moment = True
                    for q in gate_qubits:
                        p0 = self.p0_map.get(tuple([q]))
                        p1 = self.p1_map.get(tuple([q]))
                        if p0 is None:
                            raise ValueError("No p0 value for qubit {}".format(q))
                        if p1 is None:
                            raise ValueError("No p1 value for qubit {}".format(q))
                        if not np.isclose(p0, 0) or not np.isclose(p1, 0):
                            new_moment.append(channels.asymmetric_bitflip(p0=p0, p1=p1).on(q))

                elif n_qubit == 1:
                    fid_1q = self.fid_1q_map.get(tuple(gate_qubits))
                    if fid_1q is not None and fid_1q > 0:
                        new_moment.append(cirq.depolarize(p=fid_1q).on(gate_qubits[0]))
                        touched.update(gate_qubits)
                elif n_qubit == 2:
                    sorted_2q = tuple(sorted(gate_qubits, key=lambda q: q.row))
                    fid_2q_tot = self.fid_2q_tot_map.get(sorted_2q)
                    # if fid_2q_tot is None:
                    #     fid_2q_tot = self.fid_2q_tot_map.get(tuple(reversed(gate_qubits)))
                    fid_2q_pure = self.fid_2q_pure_map.get(sorted_2q)
                    # if fid_2q_pure is None:
                    #     fid_2q_pure = self.fid_2q_pure_map.get(tuple(reversed(gate_qubits)))
                    if fid_2q_pure is None or fid_2q_tot is None:
                        raise ValueError("Did not receive fid_2q value for qubits {}".format(sorted_2q))

                    # incoherent noise is a symmetric generalized 2q depol
                    # NOTE: `pure` here refers to error in the purity i.e.
                    # incoherent noise
                    if fid_2q_pure > 0:
                        # incoherent noise -> nonunitary
                        new_moment.append(
                            channels.symmetric_depolarize(p=fid_2q_pure).on(*gate_qubits))
                        touched.update(gate_qubits)

                    # coherent noise is a random 2-qubit rotation
                    fid_2q_unitary = fid_2q_tot - fid_2q_pure
                    if fid_2q_unitary > 0:
                        new_moment.append(
                            channels.random_2q_pauli_exp(p=fid_2q_unitary, qubits=gate_qubits))
                        touched.update(gate_qubits)

            new_op_tree = [new_moment]
            # Insert decohering channels according to gate durations
            # DONT apply these to readouts, as they were already touched by ABF
            if not last_moment and not np.isclose(moment_duration, 0.):
                # Apply time-decay decoherence for untouched qubits this moment
                time_decay_op_tree = []
                for q in set(circuit.all_qubits()).difference(touched):
                    T1, T2 = self.T1_map.get(tuple([q])), self.T2_map.get(tuple([q]))
                    p_decay, p_dephasing = channels.decay_prob(T1, T2, moment_duration)
                    # separate into coherent and incoherent noise
                    if not np.isclose(T1, 0):
                        time_decay_op_tree.append(cirq.amplitude_damp(gamma=p_decay).on(q))
                    if not np.isclose(T2, 0):
                        time_decay_op_tree.append(cirq.phase_damp(gamma=p_dephasing).on(q))
                new_op_tree.append(time_decay_op_tree)
            # Original gates always represent _last_ moment, in case of measurement
            new_op_tree.append(moment)
            out.append(new_op_tree)

        # If no measurement gates were found, simulate bitflip channels as the
        # last layer of the circuit. The `last_moment` flag only raises if
        # measurement gates are seen
        if not last_moment:
            new_moment = []
            for q in circuit.all_qubits():
                p0 = self.p0_map.get(tuple([q]))
                p1 = self.p1_map.get(tuple([q]))
                if not np.isclose(p0, 0) or not np.isclose(p1, 0):
                    new_moment.append(channels.asymmetric_bitflip(p0=p0, p1=p1).on(q))
            out.append(new_moment)

        return out

