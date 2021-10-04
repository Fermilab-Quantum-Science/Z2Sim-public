import pytest
import os

import numpy as np
import matplotlib.pyplot as plt
import cirq

from simple_noise_model import SimpleNoiseModel

# Toggle this on if you want to see T2 decay curves. currently, visual
# inspection shows that the test works.
DISPLAY = False

def test_noisemodel_from_fixed_values():
    noise_model = SimpleNoiseModel.from_fixed_values(
        T1=1,
        T2=2,
        p0=3,
        p1=4,
        fid_1q=5,
        fid_2q_tot=6,
        fid_2q_pure=7)
    assert noise_model.T1_map[cirq.LineQubit(0)] == 1
    assert noise_model.T2_map["dummy string"] == 2
    assert noise_model.p0_map.get(4210) == 3
    assert noise_model.p1_map.get((cirq.LineQubit(1), cirq.LineQubit(2))) == 4
    assert noise_model.fid_1q_map.get(("a", "b", "c")) == 5
    assert noise_model.fid_2q_tot_map.get(None) == 6
    assert noise_model.fid_2q_pure_map.get(None) == 7


def test_t1_free_decay():
    """Description of experiment: Starting from the state |111>, flip to/from
    |000> an even number of times then record the probability of being in |111>.
    Use either a T1 or T2 decay model to simulate an exponential probability
    decay."""
    qubits = cirq.GridQubit.rect(1, 3)
    # To get T1 decay with _fixed_ gate times, we require gate_time << T1
    # and relatively large `n_layers`
    t1_lst = [100, 200, 300]
    n_layers = 30
    repetitions = 2000

    circuit = cirq.Circuit()
    for g in range(n_layers):
        circuit += cirq.Circuit(*[cirq.X(q) for q in qubits])

    # Noise model will apply T1 style decay for the time n_gates * T(gate)
    t1_noise_only = SimpleNoiseModel.from_fixed_values()
    qubit_keys = [tuple([q]) for q in qubits]
    t1_noise_only.T1_map = dict(zip(qubit_keys, t1_lst))
    circuit = cirq.Circuit(*circuit.with_noise(t1_noise_only).all_operations(),
                                     strategy=cirq.InsertStrategy.EARLIEST)

    initial_state = np.zeros((2 ** len(qubits),) * 2, dtype=np.complex64)
    initial_state[-1,-1] = 1
    sim = cirq.DensityMatrixSimulator()
    final_state = sim.simulate(circuit, initial_state=initial_state).final_density_matrix
    results = cirq.sample_density_matrix(final_state,
                                         range(len(qubits)),
                                         repetitions=repetitions)
    excited_state_probs = np.zeros((len(qubits)))
    for j in range(len(qubits)):
        excited_state_probs[j] = np.mean(results[:, j])

    t_tot = 0.025 * n_layers
    # No factor of 2 in probabilities, since density matrix diagonals are
    # probabilities being decayed by exp(t/T1)
    truths = np.array([np.exp(-1 * t_tot/t1) for t1 in t1_lst])
    np.testing.assert_allclose(truths, excited_state_probs, atol=0.01)


@pytest.mark.skipif(not DISPLAY, reason="Currently only a visual check.")
def test_t2_free_decay():
    """Visually inspect the effects of T2 decay.

    Starting from the state |000>, rotate to
    |+++>, evolve for some amount of time, and flip to |111>. Record exponential
    decay in P(111).
    """
    qubits = cirq.GridQubit.rect(1, 3)
    # To get T1 decay with _fixed_ gate times, we require gate_time << T1
    # and relatively large `n_layers`
    t2_lst = [1, 2, 4]
    n_layers = 40
    n_layers_arr = np.arange(1, n_layers)
    all_times = 0.025 * n_layers_arr
    all_p1 = []
    for n_layers in n_layers_arr:
        repetitions = 5000

        circuit = cirq.Circuit(*[cirq.ry(np.pi/2)(q) for q in qubits])
        for n in range(n_layers):
            circuit += cirq.Circuit(*[cirq.I(q) for q in qubits])
        circuit += cirq.Circuit(*[cirq.ry(np.pi/2)(q) for q in qubits])
        circuit.append(cirq.measure(*qubits, key='m'))
        # Noise model will apply T2 style decay for the time n_gates * T(gate)
        t2_noise_only = SimpleNoiseModel.from_fixed_values()
        qubit_keys = [tuple([q]) for q in qubits]
        t2_noise_only.T2_map = dict(zip(qubit_keys, t2_lst))
        circuit = circuit.with_noise(t2_noise_only)

        trial = cirq.DensityMatrixSimulator().run(circuit, repetitions=repetitions)
        results = trial.measurements['m']
        excited_state_probs = np.zeros((len(qubits)))
        for j in range(len(qubits)):
            excited_state_probs[j] = np.mean(results[:, j])

        all_p1.append(excited_state_probs)

    truths = np.array([0.5 + 0.5 * np.exp(-0.5 * all_times/T2) for T2 in t2_lst]).T
    # import pdb; pdb.set_trace()
    colors = ['r', 'b', 'g']
    fig, ax = plt.subplots()
    # analytical expectations dashed lines
    for col, color in zip(np.array(truths).T, colors):
        ax.plot(all_times, col, c=color, ls='--')
    for col, color in zip(np.array(all_p1).T, colors):
        ax.plot(all_times, col, c=color)
    plt.show()

    # np.testing.assert_allclose(truths, excited_state_probs, atol=0.01)


@pytest.mark.parametrize('qubit_p0p1_lst', (
[
    [(0.3, 0.0), (0.5, 0.0), (0.1, 0.0)], # p0 only tests
    [(0.0, 0.5), (0.0, 0.2), (0.0, 0.1)], # p1 only tests
    [(0.1, 0.5), (0.6, 0.2), (0.1, 0.1)], # p0 and p1 tests
])
)
def test_p0_p1_statistics(qubit_p0p1_lst):
    """Description of experiment: Measure bitstrings from a blank circuit
    and confirm that their distributions are correct w/r to p0 or p1."""

    qubits = cirq.GridQubit.rect(1, 3)
    repetitions = 10000

    qubit_keys = [tuple([q]) for q in qubits]
    p0_map = dict(zip(qubit_keys, [tupl[0] for tupl in qubit_p0p1_lst]))
    p1_map = dict(zip(qubit_keys, [tupl[1] for tupl in qubit_p0p1_lst]))

    p0p2_noise = SimpleNoiseModel.from_fixed_values()
    p0p2_noise.p0_map = p0_map
    p0p2_noise.p1_map = p1_map
    # 50/50 noisy circuit with measurement errors
    circuit = cirq.Circuit(*[cirq.H(q) for q in qubits])
    circuit += cirq.Circuit(*[cirq.measure(q, key=i) for i, q in enumerate(qubits)])
    circuit = cirq.Circuit(*circuit.with_noise(p0p2_noise).all_operations(),
                                     strategy=cirq.InsertStrategy.EARLIEST)
    measurements = cirq.DensityMatrixSimulator().run(circuit, repetitions=repetitions).measurements
    ordered_measurements = np.asarray([np.asarray(measurements.get(i)) for i, _ in enumerate(qubits)])

    print(circuit)
    for i, qubit_res in enumerate(ordered_measurements):
        p0, p1 = qubit_p0p1_lst[i]
        n0 = len(qubit_res[np.where(qubit_res == 0)])
        n1 = len(qubit_res[np.where(qubit_res == 1)])

        n0_truth = int(repetitions * (1 - p0 + p1)/ 2)
        n1_truth = int(repetitions * (1 - p1 + p0)/ 2)
        err0 = (n0_truth - n0) / n0_truth
        err1 = (n1_truth - n1) / n1_truth
        assert err0 < 0.05 and err1 < 0.05


def test_xeb_no_double_touch():
    q0, q1, q2 = cirq.GridQubit.rect(3, 1)
    noise_model = SimpleNoiseModel.from_fixed_values(T1=10.)
    noise_model.fid_2q_tot_map = {(q0, q1): .5}
    noise_model.fid_2q_pure_map = {(q0, q1): .3}
    circuit = cirq.Circuit(cirq.FSimGate(np.pi/2, np.pi/6)(q0, q1))
    circuit.append(cirq.X(q2))
    noisy_circuit = circuit.with_noise(noise_model)
    # no extra moments due to T1
    # depth 3 = FSIM gate, then unitary noise, then decoherent 2q noise
    assert(len(noisy_circuit) == 3)

    noise_model = SimpleNoiseModel.from_fixed_values(T1=10.)
    noise_model.fid_1q_map = {(q0): .5}
    circuit = cirq.Circuit(*[cirq.X(q) for q in (q0, q1, q2)])
    noisy_circuit = circuit.with_noise(noise_model)
    # depth 2 = gate, then gate infidelity noise
    assert(len(noisy_circuit) == 2)


@pytest.mark.skipif(not DISPLAY, reason="Currently only a visual check.")
@pytest.mark.parametrize('p', [0.1, 0.2])
def test_1q_rb(p):
    """Implement randomized benchmarking to determine 1q rb parameters work.

    # WARNING: this takes a hecking long time to run.
    """
    q = cirq.GridQubit(0, 0)
    rb1q_gate_noise = SimpleNoiseModel.from_fixed_values(fid_1q=p)
    # Track the gate that will undo all of the applied gates
    c_r = cirq.PauliString(cirq.I(q))
    small_c1 = [cirq.X, cirq.Y, cirq.Z]
    repetitions = 1000
    # iterate over clifford sequence depths
    p_0_arr = []
    truth_arr = []
    m_max = 20
    for m in range(1,m_max):
        cliffs = np.random.choice(small_c1, size=m)
        c_r = cirq.PauliString(cliffs[-1](q))
        for c in list(reversed(cliffs))[1:]:
            c_r *= cirq.PauliString(c(q))

        circuit = cirq.Circuit(*(c(q) for c in cliffs))
        # conversion to gate
        try:
            c_r = list(c_r.values())[0].on(*c_r.qubits)
        except IndexError:
            c_r = cirq.I(q)
        circuit += c_r
        circuit += cirq.Circuit(cirq.measure(q, key='q'))
        # results = cirq.Simulator().run(program=circuit, repetitions=repetitions).measurements.get('q')
        # p_0 = len(results[results == 0]) / repetitions
        # # sanity check for noiseless C1 sequence
        # assert p_0 == 1
        noisy_circuit = circuit.with_noise(rb1q_gate_noise)
        results = cirq.Simulator().run(program=noisy_circuit, repetitions=repetitions).measurements.get('q')
        p_0_rb = len(results[results == 0]) / repetitions

        p_0_arr.append(p_0_rb)
        # An exponential converging to 1/2
        truth_arr.append(0.5 * (1-p) ** m + 0.5)

    plt.plot(range(1,m_max), p_0_arr, range(1,m_max), truth_arr)
    plt.show()


@pytest.mark.skip(reason="2q gate depolarization is not a well-studied model")
def test_2q_fid():
    # TODO: I don't know how to test these.
    pass
