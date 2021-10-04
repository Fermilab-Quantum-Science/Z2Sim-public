import numpy as np

import cirq


# TODO: Depending on the status of their serializer and what kinds of
# symbolic gates we need to implement, the following workaround may still
# be necessary...
def serializable_rz(x):
    return cirq.ZPowGate(exponent=x)


def serializable_ry(x):
    return cirq.YPowGate(exponent=x)



def CNOT_using_sqrtISWAP(q0, q1):
    """Return an operation executing CNOT(q0, q1) using the sqrtISWAP entangler.

    TODO: Some rx, rz combinations could be collapsed into PhasedXZ's but
        this is unlikely to be necessary.
        
    Args:
        q0: Control qubit
        q1: Target qubit

    Returns:
        Generator for a cirq.Operation implementing CNOT(q0, q1)
    """

    yield cirq.decompose(cirq.H(q1))
    yield cirq.rx(np.pi*0.5).on(q0)
    yield cirq.rx(np.pi*0.5).on(q1)
    yield cirq.ISWAP.on(q0, q1) ** 0.5
    yield cirq.rx(np.pi*-1.0).on(q0)
    yield cirq.Z(q0)
    yield cirq.ISWAP.on(q0, q1) ** 0.5
    yield cirq.Z(q0)
    yield cirq.rx(np.pi*-0.5).on(q1)
    yield cirq.rx(np.pi*0.5).on(q0)
    yield cirq.rz(np.pi*0.5).on(q0)
    yield cirq.rz(np.pi*0.5).on(q1)
    yield cirq.decompose(cirq.H(q1))


def CNOT_using_sycamore(q0, q1):
    """Return an operation executing CNOT(q0, q1) using the sycamore entangler.

    Args:
        q0: Control qubit
        q1: Target qubit

    Returns:
        Generator for a cirq.Operation implementing CNOT(q0, q1)
    """

    sycamore_gate = cirq.FSimGate(np.pi/2, np.pi/6).on(q0, q1)
    yield cirq.decompose(cirq.H(q1))
    yield cirq.rx(np.pi*1.4136540654256826).on(q1)
    yield sycamore_gate
    yield cirq.rz(np.pi/12).on(q0)
    yield cirq.rz(np.pi/12).on(q1)
    yield cirq.rx(np.pi*-0.4771266984986656).on(q0)
    yield cirq.Z(q0)
    yield sycamore_gate
    yield cirq.rz(np.pi/12).on(q0)
    yield cirq.rz(np.pi/12).on(q1)
    yield cirq.Z(q0)
    yield cirq.rx(np.pi*-1.4136540654256826).on(q1)
    yield cirq.rz(np.pi*0.5).on(q0), cirq.rz(np.pi*0.5).on(q1)
    yield cirq.decompose(cirq.H(q1))


def ZZ_using_param_ISWAP(q0, q1, angle):
    """
    Return a generator that executes e^{i angle ZZ} using the parameterized ISWAP gate 
    Args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementings e^{i angle ZZ}
    """
    iswap = cirq.ISWAP.on(q0, q1) ** (angle / np.pi * 2)
    yield cirq.decompose(cirq.H(q0))
    yield cirq.decompose(cirq.H(q1))
    yield cirq.X(q0)
    yield iswap
    yield cirq.X(q0)
    yield iswap
    yield cirq.decompose(cirq.H(q0))
    yield cirq.decompose(cirq.H(q1))


def ZZ_using_sqrtISWAP(q0, q1, angle):
    """
    Return a generator that executes e^{i angle ZZ} using the sqrtISWAP gate 
    Args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementings e^{i angle ZZ}
    """
    yield CNOT_using_sqrtISWAP(q0, q1)
    yield cirq.rz(-angle * 2).on(q1)
    yield CNOT_using_sqrtISWAP(q0, q1)


def XX_using_sqrtISWAP(q0, q1, angle):
    """
    Return a generator that executes e^{i angle ZZ} using the sqrtISWAP gate 
    Args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementings e^{i angle ZZ}
    """
    yield CNOT_using_sqrtISWAP(q0, q1)
    yield cirq.rx(-angle * 2).on(q0)
    yield CNOT_using_sqrtISWAP(q0, q1)

    
def XX_using_param_ISWAP(q0, q1, angle):
    """
    Return a generator that executes e^{i angle XX} using the
    parameterized ISWAP gate 
    Args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementings e^{i angle XX}
    """
    iswap = cirq.ISWAP.on(q0, q1) ** (angle / np.pi * 2)
    yield cirq.X(q0)
    yield iswap
    yield cirq.X(q0)
    yield iswap
    

def ZZ_using_param_sycamore(q0, q1, angle):
    """
    Return a generator that executes e^{i angle ZZ} using the
    parameterized sycamore gate 
    Args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementings e^{i angle ZZ}
    """
    sycamore_gate = cirq.FSimGate(0, -angle * 4).on(q0, q1)
    yield sycamore_gate
    yield cirq.Z(q0) ** (-angle * 2 / np.pi)
    yield cirq.Z(q1) ** (-angle * 2 / np.pi)


def ZZ_using_sycamore(q0, q1, angle):
    """
    Return a generator that executes e^{i angle ZZ} using the 
    default sycamore gate 
    Args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementings e^{i angle ZZ}
    """
    yield CNOT_using_sycamore(q0, q1)
    yield cirq.Z(q1) ** (-angle * 2 / np.pi)
    yield CNOT_using_sycamore(q0, q1)


def XX_using_param_sycamore(q0, q1, angle):
    """
    Return a generator that executes e^{i angle XX} using the 
    parameterized sycamore gate 
    Args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementings e^{i angle XX}
    """
    sycamore_gate = cirq.FSimGate(0, -angle * 4).on(q0, q1)
    yield cirq.decompose(cirq.H(q0))
    yield cirq.decompose(cirq.H(q1))
    yield sycamore_gate
    yield cirq.Z(q0) ** (-angle * 2 / np.pi)
    yield cirq.Z(q1) ** (-angle * 2 / np.pi)
    yield cirq.decompose(cirq.H(q0))
    yield cirq.decompose(cirq.H(q1))


def XX_using_sycamore(q0, q1, angle):
    """
    Return a generator that executes e^{i angle XX} using the default
    sycamore gate 
    Args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementings e^{i angle XX}
    """
    yield CNOT_using_sycamore(q0, q1)
    yield cirq.rx(-angle * 2).on(q0)
    yield CNOT_using_sycamore(q0, q1)


def XX_using_CNOT(q0, q1, angle):
    """
    Return a generator that executes e^{i angle XX} using a CNOT gate
    args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementing e^{i angle XX}
    """
    yield cirq.CNOT(q0, q1)
    yield cirq.X(q0) ** (-angle * 2 / np.pi)
    yield cirq.CNOT(q0, q1)

    
def ZZ_using_CNOT(q0, q1, angle):
    """
    Return a generator that executes e^{i angle XX} using a CNOT gate
    args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementing e^{i angle XX}
    """
    yield cirq.CNOT(q0, q1)
    yield cirq.rz(-angle * 2).on(q1)
    yield cirq.CNOT(q0, q1)
    

    
def XX_using_CZ(q0, q1, angle):
    """
    Return a generator that executes e^{i angle XX} using a CZ gate
    args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementing e^{i angle XX}
    """
    yield cirq.decompose(cirq.H(q1))
    yield cirq.CZ(q0, q1)
    yield cirq.rx(-angle * 2).on(q0)
    yield cirq.CZ(q0, q1)
    yield cirq.decompose(cirq.H(q1))


def ZZ_using_CZ(q0, q1, angle):
    """
    Return a generator that executes e^{i angle ZZ} using a CNOT gate
    args:
        q0: qubit 0
        q1: qubit 1
    
    Returns:
        Generator for a cirq.Operation implementing e^{i angle ZZ}
    """
    yield cirq.decompose(cirq.H(q1))
    yield cirq.CZ(q0, q1)
    yield cirq.rx(-angle * 2).on(q1)
    yield cirq.CZ(q0, q1)
    yield cirq.decompose(cirq.H(q1))
    
