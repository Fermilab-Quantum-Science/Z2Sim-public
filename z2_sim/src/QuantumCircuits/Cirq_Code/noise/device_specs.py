"""Hardware-dependent parameters for tuning a noise model."""

import numpy as np
import cirq


def get_gate_duration(op: cirq.ops.gate_operation.GateOperation,
                      atol: float = 1e-4, verbose=True) -> float:
    """Look up the duration (microseconds) for a gate instance.

    These are based on my most recent knowledge of the Sycamore "Rainbow" chip,
    from sometime in early 2020.

    To be an instance, `op` should have a valid exponent and be initialized
    on a qubit(s). If a two qubit gate instance is passed, it must be a valid
    Sycamore hardware gate (ISwap ** 0.5 or FSimGate ** x).

    Args:
        A cirq GateOperation to look up execution time for

    Returns:
        gate execution time in microseconds.
    """

    if not hasattr(op, "gate"):
        raise TypeError("Expected an initialized gate with a `gate` member. "
                        "Instead got {} with type {}."
                        " ".format(op, type(op)))
    gate = op.gate

    if isinstance(gate, cirq.ops.identity.IdentityGate):
        # For noise model tests
        return 0.025

    if isinstance(gate, cirq.ops.pauli_gates._PauliZ):
        return 0.0  # Z gate is just a reference fram switch with post-processing

    if isinstance(gate, (cirq.SingleQubitGate, cirq.IdentityGate)):
        return 0.025

    if isinstance(gate, cirq.ISwapPowGate) and abs(gate.exponent) == 0.5:
        return 0.032  # FIXME: DOES EXPONENT = -0.5 COUNT?

    if isinstance(gate, cirq.FSimGate) and (np.isclose(gate.theta, np.pi / 2)
                                            and np.isclose(gate.phi, np.pi / 6)):
        return 0.012

    if isinstance(gate, (cirq.CXPowGate, cirq.CZPowGate)):
        # Here I'll assume a fsim decomposition
        return 5 * 0.025 + 2 * 0.012

    if isinstance(gate, cirq.MeasurementGate):
        return 4.0

    # Patch for gates outside our hardware gateset:
    if verbose:
        print("Recieved gate {}, which does not have a valid "
              "duration for hardware runs. Using default times".format(op))
    if gate.num_qubits() == 1:
        return 0.015

    if gate.num_qubits() == 2:
        return 0.03

    # Should be out of reach
    raise ValueError("Recieved gate {}, which does not have a valid "
                     "duration for hardware runs.".format(op))