import cirq


def make_all_plaq_observables(source: cirq.GridQubit, ancilla: cirq.GridQubit, circuit: cirq.Circuit):
    """Generate a list of observables like

        <X_i (X + iY)_ancilla >

    for every location `i` in the grid.

    Args:
        source, ancilla: "special" qubits for this correlator circuit
        circuit: Correlator simulation circuit
    """
    # Observable is \prod_{k\in grid} X_k (X + iY)_ancilla
    gridqubits = list(set(circuit.all_qubits()))
    gridqubits.remove(ancilla)
    observables = []
    for src in gridqubits:
        observables.append(cirq.X(src) * (cirq.X(ancilla) + 1j * cirq.Y(ancilla)) )
    return observables