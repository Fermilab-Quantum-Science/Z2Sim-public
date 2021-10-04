"""v1 - Try just replicating some noiseless sweep data taking runs."""
from datetime import datetime
import time

import qsimcirq
import cirq
import numpy as np

from z2_sim.src.QuantumCircuits.Cirq_Code.Z2GaugeCirq import make_trotter_circuit
from z2_sim.src.QuantumCircuits.Cirq_Code import util


grid_sizes = [3, 4, 5]
N_TROTTER_STEPS = 20
EFIELD = 10
MFIELD = 1 / EFIELD
DT = 0.05
obc = True

N_THREADS = 8
FUSION = 4   # previously optimized
TRAJECTORIES = 10
sim_types = [
    {'r': TRAJECTORIES, 'f': FUSION, 'g': True}, # GPU
    {'r': TRAJECTORIES, 'f': FUSION, 't': N_THREADS}, # CPU
    ]

timestamp = datetime.today().strftime('%Y%m%d')
time_results = np.zeros((2, len(grid_sizes)))

for i, qsim_options in enumerate(sim_types):
    if i == 0:
        print("GPU:")
    else:
        print("CPU")
    for j, n in enumerate(grid_sizes):
        print("\t grid size: {}".format(n))
        qsim_simulator = qsimcirq.QSimSimulator(qsim_options=qsim_options)
        circuit, (source, ancilla) = make_trotter_circuit(n, n, n_timesteps=N_TROTTER_STEPS, Jcoup=EFIELD, Gamma=MFIELD, dt=DT, obc=obc, decompbasis='MS')
        noisy_circuit = circuit.with_noise(cirq.phase_damp(0.01))
        observables = util.make_all_plaq_observables(source, ancilla, circuit)
        start = time.time()
        expectation = qsim_simulator.simulate_expectation_values_sweep(noisy_circuit, observables=observables, params=None)[0]
        end = time.time()

        time_results[i, j] = end - start
        print("\t", end - start)

np.save(f"results_gpu_bench_{timestamp}.npy", time_results)