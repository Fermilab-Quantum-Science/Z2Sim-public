"""FIXME: update for new module structure if ever necessary."""
import cirq
from numpy.random import rand
from numpy import pi, cos, sin, identity, array

from context import decompositions


DEVICE = cirq.google.devices.Sycamore23
# Here we don't care about the choice of qubits, only that the operations
# are valid.
TEST_QUBITS = (cirq.GridQubit(4, 1), cirq.GridQubit(4, 2))


def test_cnot_from_sqrtiswap():
    q0, q1 = cirq.LineQubit.range(2)
    unitary = cirq.unitary(
        cirq.Circuit(*decompositions.CNOT_using_sqrtISWAP(q0, q1)))
    truth = cirq.unitary(cirq.CNOT(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_cnot_from_sqrtiswap_compatibility():
    circuit = cirq.Circuit(*decompositions.CNOT_using_sqrtISWAP(*TEST_QUBITS))
    DEVICE.validate_circuit(circuit)


def test_cnot_from_sycamore():
    q0, q1 = cirq.LineQubit.range(2)
    unitary = cirq.unitary(
        cirq.Circuit(*decompositions.CNOT_using_sycamore(q0, q1)))
    truth = cirq.unitary(cirq.CNOT(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_cnot_from_sycamore_compatibility():
    circuit = cirq.Circuit(*decompositions.CNOT_using_sycamore(*TEST_QUBITS))
    DEVICE.validate_circuit(circuit)


def test_xx_using_param_iswap():
    q0, q1 = cirq.LineQubit.range(2)
    angle = rand()
    unitary = cirq.unitary(cirq.Circuit(*decompositions.XX_using_param_ISWAP(q0, q1, angle)))
    truth = cirq.unitary(cirq.XXPowGate(exponent=-angle / pi * 2).on(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_xx_using_sqrt_iswap():
    q0, q1 = cirq.LineQubit.range(2)
    angle = rand()
    unitary = cirq.unitary(cirq.Circuit(*decompositions.XX_using_sqrtISWAP(q0, q1, angle)))
    truth = cirq.unitary(cirq.XXPowGate(exponent=-angle / pi * 2).on(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_zz_using_param_iswap():
    q0, q1 = cirq.LineQubit.range(2)
    angle = rand()
    unitary = cirq.unitary(cirq.Circuit(*decompositions.ZZ_using_param_ISWAP(q0, q1, angle)))
    truth = cirq.unitary(cirq.ZZPowGate(exponent=-angle / pi * 2).on(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_zz_using_sqrt_iswap():
    q0, q1 = cirq.LineQubit.range(2)
    angle = rand()
    unitary = cirq.unitary(cirq.Circuit(*decompositions.ZZ_using_sqrtISWAP(q0, q1, angle)))
    truth = cirq.unitary(cirq.ZZPowGate(exponent=-angle / pi * 2).on(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_xx_using_param_sycamore():
    q0, q1 = cirq.LineQubit.range(2)
    angle = rand()
    unitary = cirq.unitary(cirq.Circuit(*decompositions.XX_using_param_sycamore(q0, q1, angle)))
    truth = cirq.unitary(cirq.XXPowGate(exponent=-angle / pi * 2).on(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_xx_using_sycamore():
    q0, q1 = cirq.LineQubit.range(2)
    angle = rand()
    unitary = cirq.unitary(cirq.Circuit(*decompositions.XX_using_sycamore(q0, q1, angle)))
    truth = cirq.unitary(cirq.XXPowGate(exponent=-angle / pi * 2).on(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_zz_using_param_sycamore():
    q0, q1 = cirq.LineQubit.range(2)
    angle = rand()
    unitary = cirq.unitary(cirq.Circuit(*decompositions.ZZ_using_param_sycamore(q0, q1, angle)))
    truth = cirq.unitary(cirq.ZZPowGate(exponent=-angle / pi * 2).on(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_zz_using_sycamore():
    q0, q1 = cirq.LineQubit.range(2)
    angle = rand()
    unitary = cirq.unitary(cirq.Circuit(*decompositions.ZZ_using_sycamore(q0, q1, angle)))
    truth = cirq.unitary(cirq.ZZPowGate(exponent=-angle / pi * 2).on(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_zz_using_cnot():
    q0, q1 = cirq.LineQubit.range(2)
    angle = rand()
    unitary = cirq.unitary(cirq.Circuit(*decompositions.ZZ_using_CNOT(q0, q1, angle)))
    truth = cirq.unitary(cirq.ZZPowGate(exponent=-angle / pi * 2).on(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)


def test_xx_using_cnot():
    q0, q1 = cirq.LineQubit.range(2)
    angle = rand()
    unitary = cirq.unitary(cirq.Circuit(*decompositions.XX_using_CNOT(q0, q1, angle)))
    truth = cirq.unitary(cirq.XXPowGate(exponent=-angle / pi * 2).on(q0, q1))
    cirq.testing.assert_allclose_up_to_global_phase(unitary, truth, atol=1e-6)

