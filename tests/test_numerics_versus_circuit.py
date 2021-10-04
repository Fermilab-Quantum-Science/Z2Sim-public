"""test_numerics_versus_circuit.py

This is a test for consistency between Cirq circuits and numpy/scipy.sparse
based implementation of smaller grids.


TODO:
    - should replace the current trotter simulation scheme with andy's
        intermediate state simulator.
"""
import numpy as np
import matplotlib.pyplot as plt
import cirq
import qsimcirq

from z2_sim.src.NumericsPython.Z2DualCorrelation import z2_dual_correlater_numeric
from z2_sim.src.QuantumCircuits.Cirq_Code.Z2GaugeCirq import make_trotter_circuit




NGRID = 3
JCOUP = 10
GAMMA = 0.1
DT = 1
N_TROTTER_STEPS = 40

# Do the truth numeric simulation
truth = z2_dual_correlater_numeric(NGRID, JCOUP, GAMMA, DT, N_TROTTER_STEPS, False, False)


# Initialize the simulator
N_REPETITIONS = 1 # Not relevant for noiseless
N_THREADS = 32
FUSION = 4 # previously optimized
qsim_simulator = qsimcirq.QSimSimulator(qsim_options={'r': N_REPETITIONS, 't': N_THREADS, 'f': FUSION})

results = np.zeros((NGRID, NGRID,N_TROTTER_STEPS + 1), dtype=complex)


for ell, n_trotter_steps in enumerate(range(1, N_TROTTER_STEPS, 1)):
    circuit, (source, ancilla) = make_trotter_circuit(NGRID, NGRID,
                                                      n_timesteps=n_trotter_steps,
                                                      Jcoup=JCOUP,
                                                      Gamma=GAMMA,
                                                      dt=DT, decompbasis='CNOT')

    qubits = list(circuit.all_qubits())
    nqubits = len(qubits)

    # Observable is \prod_{k\in grid} X_k (X + iY)_ancilla
    gridqubits = set(circuit.all_qubits())
    gridqubits.remove(ancilla)
    observables = []
    for src in gridqubits:
        observables.append(cirq.X(src) * (cirq.X(ancilla) + 1j * cirq.Y(ancilla)) )
    expectations = qsim_simulator.simulate_expectation_values_sweep(
        circuit,
        observables=observables,
        params=None,
        permit_terminal_measurements=True)[0] # since sweep is sz 1

    # Now assign these expectation values to the results
    # CAREFUL: nx = number of rows in this case.
    for idx, q in enumerate(gridqubits):
        results[q.row, q.col, ell+1] = expectations[idx]

plt.plot(range(N_TROTTER_STEPS+1), results.mean(axis=(0,1)), label="cirq code" )
plt.plot(range(N_TROTTER_STEPS+1), truth.mean(axis=(0,1)), label="numpy numerics")
plt.xlabel("trotter step")
plt.ylabel("Average <X_i (X + iY)_anc>")
plt.legend()
plt.show()