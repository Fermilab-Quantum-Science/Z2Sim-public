"""production.py - production code for computing observables."""
import time

from typing import Optional, Callable, Sequence

import numpy as np
import cirq
import qsimcirq

from z2_sim.src.QuantumCircuits.Cirq_Code.Z2GaugeCirq import (
    make_trotter_circuit,
    make_trotter_circuit_ancillafree,
    make_ancillafree_inputs
)
from z2_sim.src.QuantumCircuits.Cirq_Code import util


def compute_obs_with_intermediate_state_vector(
    n: int,
    trotter_steps: int,
    jcoup: float,
    dt: float,
    all_observables: Optional[bool] = True,
    qsim_options: Optional[dict] = dict(t=32, f=4),
    decompbasis: Optional[str] = None,
    obc: Optional[bool] = True,
    noise_models: Optional[Sequence[cirq.devices.NoiseModel]] = None,
    verbose: Optional[bool] = False,
) -> np.ndarray:
    """
    Compute the time-dependent expectation value

        f_i = <0|U*(t) X_i U(t) X_s|0>

    for the Z2 trotter simulation.

    TODO:
        - units for J, dt?
        - Add option to combine single qubit gates and Ising exp(isZZ) style gates?
            Need to consult with Martin about the representation for these gates
            as they will likely not be supported with out-of-the-box cirq.

    Args:
        n: dimension of grid. Total number of qubits is `n ** 2 + 1`
        trotter_steps: Total number of trotter steps to simulate, with observable
            sampled after each step.
        jcoup: Transverse/E-field parameter.
        dt: Timestep size.
        all_observables: If True, compute the expectation value `f_i` for every
            i=1, ..., n**2. Otherwise, compute a single expectation value at
            each step using a source qubit near the center of the grid.
        decompbasis: Decomposition basis for the trotter circuits. If `None`,
            this will use `exp(i * m * ZZ)` Ising style gates natively
        obc: Flag setting open boundary conditions on the trotter simulation.
        verbose: Print statements to stdout.

    Returns:
        If `all_observables` is True, return an array with shape

            `(trotter_steps, n, n)`

        Otherwise, return an array with shape `(trotter_steps,)`. The first
        dimension of these outputs contains the observable values corresponding
        to each timestep from `i=1,...,trotter_steps`.
    """

    vprint = lambda *args: print(*args) if verbose else None
    if decompbasis is None:
        decompbasis = 'MS'

    if verbose:
        # store times for computing expectation values and simulating
        # first axis: each trotter step
        # second axis: 0th element is computation, 1st element is observable time
        times_out = np.zeros((trotter_steps, 2))

    # Generate the first timestep that also prepares the correlator.
    initial_circuit, (source, ancilla) = make_trotter_circuit(n, n, n_timesteps=1, Jcoup=jcoup, Gamma=1/jcoup, dt=dt, initial=True, obc=obc, decompbasis=decompbasis)
    # Generate the unitaries applying all remaining trottersteps after the first.
    stepper_circuit, (source, ancilla) = make_trotter_circuit(n, n, n_timesteps=1, Jcoup=jcoup, Gamma=1/jcoup, dt=dt, initial=False, obc=obc, decompbasis=decompbasis)
    # In the case of noisy simulations we consider the contribution of noise
    # in the state preparation circuit to be negligible compared to the effects
    # due to noisy trotter steps.
    # WARNING: This will "stack" noise models; make sure that the noise models
    # all tag their noisy operations with the appropriate tag!!
    if noise_models:
        for noise_model in noise_models:
            stepper_circuit = stepper_circuit.with_noise(noise_model)
    qsim_simulator = qsimcirq.QSimSimulator(qsim_options=qsim_options)

    vprint("Number of qubits:", len(stepper_circuit.all_qubits()))
    vprint("\tJ={}, dt={}".format(jcoup, dt))
    observables = [cirq.X(source) * cirq.X(ancilla) + 1j * cirq.X(source) * cirq.Y(ancilla)]
    # Mirror observables over the diagonal
    if all_observables:
        triu_idx = []
        for i in range(n):
            for j in range(n):
                if j < i:
                    continue
                triu_idx.append((i, j))

        observables = [cirq.X(cirq.GridQubit(*idx)) * cirq.X(ancilla) + 1j * cirq.X(cirq.GridQubit(*idx)) * cirq.Y(ancilla) for idx in triu_idx]

    # Here we save intermediate wavefunctions and compute observables by
    # invoking qsim's function on an empty circuit with the saved state as input
    empty_c = cirq.Circuit()
    qubits = stepper_circuit.all_qubits()
    empty_c += cirq.IdentityGate(len(qubits)).on(*qubits)

    out = np.zeros((trotter_steps), dtype=np.complex64)
    if all_observables:
        out = np.zeros((trotter_steps, n, n), dtype=np.complex64)

    state = None
    for k in range(trotter_steps):
        # The first trotter step has a distinct generator for an initial state
        circuit = stepper_circuit
        if k == 0:
            circuit = initial_circuit
        t0 = time.time()

        result = qsim_simulator.simulate_sweep(
            circuit, initial_state=state, params=None)[0]

        dt = time.time() - t0
        if verbose:
            print("computation dt: ", dt)
            times_out[k, 0] = dt
        t0 = time.time()

        state = result.final_state_vector
        expectation_batch = qsim_simulator.simulate_expectation_values_sweep(
            empty_c, observables=observables, initial_state=state,
            params=None)[0]
        dt = time.time() - t0
        if verbose:
            print("observable dt: ", dt)
            times_out[k, 1] = dt
        if all_observables:
            for idx, (row, col) in enumerate(triu_idx):
                out[k, row, col] = expectation_batch[idx]
                out[k, col, row] = expectation_batch[idx]
        else:
            out[k] = expectation_batch[0]
    if verbose:
        return out, times_out

    return out


def compute_ancillafree_obs_with_intermediate_state_vector(
    n: int,
    trotter_steps: int,
    jcoup: float,
    dt: float,
    all_observables: Optional[bool] = True,
    qsim_options: Optional[dict] = dict(t=32, f=4),
    decompbasis: Optional[str] = None,
    obc: Optional[bool] = True,
    noise_models: Optional[Sequence[cirq.devices.NoiseModel]] = None,
    verbose: Optional[bool] = False,
) -> np.ndarray:
    """
    Compute the time-dependent expectation value

        f_i = <0|U*(t) X_i U(t) X_s|0>

    for the Z2 trotter simulation without using any ancilla.

    For documentation see `compute_obs_with_intermediate_state_vector`
    """

    vprint = lambda *args: print(*args) if verbose else None
    if decompbasis is None:
        decompbasis = 'MS'

    if verbose:
        # store times for computing expectation values and simulating
        # first axis: each of the 4 pseudo-observables
        # second axi: each trotter step
        # third axis: 0th element is computation, 1st element is observable time
        times_out = np.zeros((4, trotter_steps, 2))

    # Generate the unitaries applying all remaining trottersteps after the first.
    stepper_circuit, source = make_trotter_circuit_ancillafree(n, n, n_timesteps=1, Jcoup=jcoup, Gamma=1/jcoup, dt=dt, obc=obc, decompbasis=decompbasis)
    # In the case of noisy simulations we consider the contribution of noise
    # in the state preparation circuit to be negligible compared to the effects
    # due to noisy trotter steps.
    # WARNING: This will "stack" noise models; make sure that the noise models
    # all tag their noisy operations with the appropriate tag!!
    if noise_models:
        for noise_model in noise_models:
            stepper_circuit = stepper_circuit.with_noise(noise_model)
    qubits = list(stepper_circuit.all_qubits())

    # Construct a list of four circuits whose outputs will combine to yield
    # the desired observable
    input_circuits = make_ancillafree_inputs(source, qubits)

    qsim_simulator = qsimcirq.QSimSimulator(qsim_options=qsim_options)

    vprint("Number of qubits:", len(stepper_circuit.all_qubits()))
    vprint("\tJ={}, dt={}".format(jcoup, dt))
    observables = [cirq.X(source)]
    # Mirror observables over the diagonal
    if all_observables:
        triu_idx = []
        for i in range(n):
            for j in range(n):
                if j < i:
                    continue
                triu_idx.append((i, j))

        observables = [cirq.X(cirq.GridQubit(*idx)) for idx in triu_idx]

    # Here we save intermediate wavefunctions and compute observables by
    # invoking qsim's function on an empty circuit with the saved state as input
    empty_c = cirq.Circuit(cirq.I.on_each(qubits))

    # We will initially compute four separate observables, and then sum along
    # the last axis to get the desired observable
    out_raw = np.zeros((4, trotter_steps), dtype=np.complex64)
    if all_observables:
        out_raw = np.zeros((4, trotter_steps, n, n), dtype=np.complex64)
    for j in range(4):
        # Select the appropriate state-preparation circuit
        initial_circuit = input_circuits[j]
        state = None
        for k in range(trotter_steps):
            # The first trotter step has a distinct generator for an initial state
            circuit = stepper_circuit
            if k == 0:
                circuit = initial_circuit + stepper_circuit

            t0 = time.time()
            result = qsim_simulator.simulate_sweep(
                circuit, initial_state=state, params=None)[0]
            dt = time.time() - t0
            if verbose:
                print("computation dt: ", dt)
                times_out[j, k, 0] = dt

            t0 = time.time()
            state = result.final_state_vector
            expectation_batch = qsim_simulator.simulate_expectation_values_sweep(
                empty_c, observables=observables, initial_state=state,
                params=None)[0]
            dt = time.time() - t0
            if verbose:
                print("observable dt: ", dt)
                times_out[j, k, 1] = dt
            if all_observables:

                for idx, (row, col) in enumerate(triu_idx):
                    out_raw[j, k, row, col] = expectation_batch[idx]
                    out_raw[j, k, col, row] = expectation_batch[idx]
            else:
                out_raw[j, k] = expectation_batch[0]

    # Postprocessing the ouputs of the circuit sweep follows from some algebra.
    # See the appendix of the manuscript draft.
    dim_arr = np.ones((1, out_raw.ndim), int).ravel()
    dim_arr[0] = -1
    mask = np.array([-1j, 1j, 1, -1]).reshape(dim_arr) / 2
    out = (out_raw * mask).sum(axis=0)
    if verbose:
        return out, times_out

    return out


def compute_noisy_obs(
    n: int,
    trotter_start: int,
    trotter_stop: int,
    jcoup: float,
    dt: float,
    all_observables: Optional[bool] = True,
    qsim_options: Optional[dict] = dict(t=32, f=4, r=100),
    decompbasis: Optional[str] = None,
    obc: Optional[bool] = True,
    noise_models: Optional[Sequence[cirq.devices.NoiseModel]] = None,
    verbose: Optional[bool] = False,
) -> np.ndarray:
    """
    Compute the time-dependent expectation value

        f_i = <0|U*(t) X_i U(t) X_s|0>

    for the Z2 trotter simulation.

    For documentation see `compute_obs_with_intermediate_state_vector`. Docs
    below detail the modified arguments. NOTE: for trotter_stop=y,
    trotter_start=x, the results will contain a total of `y - x` trotter steps
    simulated. Runtime scales like `poly(y - x)`.

    Args:
        trotter_start: At which trotter step to begin the simulation. The first
            entry will be the output of a circuit containing `trotter_start`
            many strotter steps (i.e. inclusive indexing).
        trotter_stop: At which trotter step to end the simulation. The final
            entry will be the output of a circuit containing `trotter_stop - 1`
            many trotter steps (i.e exclusive indexing).

    Returns:
        If `all_observables` is True, return an array with shape

            `(trotter_stop - trotter_start, n, n)`
    """

    vprint = lambda *args: print(*args) if verbose else None
    if decompbasis is None:
        decompbasis = 'MS'

    if trotter_start < 1:
        raise ValueError("`trotter_start` must be at least 1.")
    if trotter_stop <= trotter_start:
        raise ValueError("`trotter_stop` must be greater than `trotter_start`")

    # Just get a template for defining observables.
    template, (source, ancilla) = make_trotter_circuit(n, n, n_timesteps=1, Jcoup=jcoup, Gamma=1/jcoup, dt=dt, initial=True, obc=obc, decompbasis=decompbasis)

    vprint("Number of qubits:", len(template.all_qubits()))
    vprint("\tJ={}, dt={}".format(jcoup, dt))
    observables = [cirq.X(source) * cirq.X(ancilla) + 1j * cirq.X(source) * cirq.Y(ancilla)]

    # Mirror observables over the diagonal
    if all_observables:
        triu_idx = []
        for i in range(n):
            for j in range(n):
                if j < i:
                    continue
                triu_idx.append((i, j))

        observables = [cirq.X(cirq.GridQubit(*idx)) * cirq.X(ancilla) + 1j * cirq.X(cirq.GridQubit(*idx)) * cirq.Y(ancilla) for idx in triu_idx]
    #
    # if all_observables:
    #     observables = util.make_all_plaq_observables(source, ancilla, template)
    #     # Reserve the identities of the source qubits for later
    #     gridqubits = list(set(template.all_qubits()))
    #     gridqubits.remove(ancilla)

    trotter_steps = trotter_stop - trotter_start
    out = np.zeros((trotter_steps), dtype=np.complex64)
    if all_observables:
        out = np.zeros((trotter_steps, n, n), dtype=np.complex64)

    # The off-by-one is due to the initial circuit containing a trotter step
    qsim_simulator = qsimcirq.QSimSimulator(qsim_options=qsim_options)

    for k, n_timesteps in enumerate(range(trotter_start, trotter_stop, 1)):

        # Generate the unitaries applying all remaining trottersteps after the first.
        complete_circuit, (source, ancilla) = make_trotter_circuit(n, n, n_timesteps=n_timesteps, Jcoup=jcoup, Gamma=1/jcoup, dt=dt, initial=True, obc=obc, decompbasis=decompbasis)
        # WARNING: This will "stack" noise models; make sure that the noise models
        # all tag their noisy operations with the appropriate tag!!
        if noise_models:
            for noise_model in noise_models:
                complete_circuit = complete_circuit.with_noise(noise_model)


        cts = count_operations_by_nqubits(complete_circuit)
        circuit_len = len(complete_circuit)
        print("RUNNING CIRCUIT FOR ")
        print(f"\tlength={circuit_len}")
        print("\t",cts)
        print(f"\ttrottersteps={n_timesteps}")
        print("\toptions:", qsim_options)

        expectation_batch = qsim_simulator.simulate_expectation_values(
            complete_circuit,
            observables=observables,
            initial_state=None,
            )

        if all_observables:
            for idx, (row, col) in enumerate(triu_idx):
                out[k, row, col] = expectation_batch[idx]
                out[k, col, row] = expectation_batch[idx]
        else:
            out[k] = expectation_batch[0]

    return out


## DIAGNOSTIC CODE: DELETEME
from collections import Counter
def count_operations_by_nqubits(circuit):
    ops = list(circuit.all_operations())
    counts = [x.gate.num_qubits() for x in ops]
    return Counter(counts)


def run_noiseless_parameter_sweep(
    n: int,
    trotter_steps: int,
    jcoup_arr: np.ndarray,
    dt_arr: np.ndarray,
    func: Callable,
    all_observables: Optional[bool] = True,
    **kwargs,
) -> np.ndarray:
    """Simulate a sweep over a set of physical parameters.

    This will compute the observables at every timestep, for every pair of
    parameters between `jcoup_arr` and `dt_arr`.

    Args:
        n: dimension of grid. Total number of qubits is `n ** 2 + 1`
        trotter_steps: Total number of trotter steps to simulate, with observable
            sampled after each step.
        jcoup_arr: Sweep over transverse/E-field parameter.
        dt_arr: Sweep over timestep size.
        func: This is the function with which to perform the sweep. This is
            expected to have a fixed signature, see the function

                `compute_obs_with_intermediate_state_vector`

            for reference.
        all_observables: If True, compute the expectation value `f_s` for every
            s=1, ..., n**2. Otherwise, compute a single expectation value at
            each step using a source qubit near the center of the grid.
        kwargs: see `compute_obs_with_intermediate_state_vector`.

    Returns:
        If `all_observables` is True, return an array with shape

            `(len(jcoup_arr), len(dt_arr), trotter_steps, n, n)`

        Otherwise, return an array with shape
            `(len(jcoup_arr), len(dt_arr), trotter_steps)`.
    """
    n_j = len(jcoup_arr)
    n_dt = len(dt_arr)

    out = np.zeros((n_j, n_dt, trotter_steps), dtype=np.complex64)
    if all_observables:
        out = np.zeros((n_j, n_dt, trotter_steps, n, n), dtype=np.complex64)

    for j, jcoup in enumerate(jcoup_arr):
        for k, dt in enumerate(dt_arr):
            out[j,k] = func(n, trotter_steps, jcoup, dt, all_observables=all_observables, **kwargs)

    return out