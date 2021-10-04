import cirq
import numpy as np

from z2_sim.src.QuantumCircuits.Cirq_Code import decompositions


def make_trotter_circuit(nx, ny, n_timesteps, Jcoup, Gamma, dt, initial=False, obc=True, decompbasis='CNOT'):
    """Generate a trotter circuit for the Z2 sim.

    This is a thin wrapper around the `Z2DualCirq` that automates selection of
    a source and ancilla qubit for computing the expectation <XX + iXY>.

    """
    efield = Jcoup
    mfield = Gamma
    z2dc = Z2DualCirq(efield, mfield, nx, ny, obc=obc, xxbasis=False, decompbasis=decompbasis)

    qubits = cirq.GridQubit.rect(nx, ny)
    # The source should be near the center. If nx,ny are odd then this will be
    # the exact center, and if n is even this will be no more than 1 cell away
    # from the center
    source = cirq.GridQubit(nx // 2, ny // 2)
    # The central source cannot be coupled with an external ancilla with
    # nearest neighbor coupling; here we'll just assume theres a small O(n)
    # overhead to do a swap network between the two.
    ancilla = cirq.GridQubit(nx, ny)
    qubits.append(ancilla)

    stateprepcirc = cirq.Circuit(cirq.I(ancilla))
    if initial:
        stateprepcirc = z2dc.generate_initial_state_correlator(source, ancilla)
    trottercircuit = z2dc.generate_trotter(n_timesteps, dt, qubits)
    return stateprepcirc + trottercircuit, (source, ancilla)


def make_trotter_circuit_ancillafree(nx, ny, n_timesteps, Jcoup, Gamma, dt, obc=True, decompbasis='CNOT'):
    """Generate a trotter circuit for the Z2 sim with no ancilla placeholder.

    This is a thin wrapper around the `Z2DualCirq` that automates selection of
    a source qubit.
    """
    z2dc = Z2DualCirq(Jcoup, Gamma, nx, ny, obc=obc, xxbasis=False, decompbasis=decompbasis)

    qubits = cirq.GridQubit.rect(nx, ny)
    source = cirq.GridQubit(nx // 2, ny // 2)
    trottercircuit = z2dc.generate_trotter(n_timesteps, dt, qubits)

    return trottercircuit, source


def make_ancillafree_inputs(source, qubits):
    """Construct the circuit inputs for the ancilla-free simulation.

    These are ordered to be consistent with the input signature of
    `util.postprocess_ancillafree_outputs`.

    """
    not_source = [x for x in qubits if x != source]
    padding = cirq.I.on_each(not_source)
    generators = [
        [cirq.rx(-np.pi/2).on(source)], # S^im+
        [cirq.rx(np.pi/2).on(source)],  # S^im-
        [cirq.H(source)],  # S^re+
        [cirq.X.on(source), cirq.H(source)] # S^re-
    ]
    out = [cirq.Circuit(*generator, padding) for generator in generators
    ]
    return out


def ZZ_native(q0, q1, angle):
    return cirq.ZZPowGate(exponent=-angle / np.pi * 2).on(q0, q1)


def XX_native(q0, q1, angle):
    return cirq.XXPowGate(exponent=-angle / np.pi * 2).on(q0, q1)


class Z2DualCirq(object):
    '''
    this class allows you to generate circuits for the Z2Gauge Theory in
    2+1 dimensions using the dual representation using Google's Cirq

    class parameters:
        nx, ny: the x and y dimension of the state
        mfield: the spacial plaquette strength (magnetic field)
        efield: the temporal plaquette strength (electric field)
        obc: a boolean flag indicating whether or not to use open boundary
            conditions. The default is true
        xxbasis: a boolean flag indicating whether to use the xxbasis for the
                nearest neighbor coupling or the zzbasis
        decompbasis: a string indicating which encoding to use for the two
                    qubit interaction valid choices:
                - 'MS': MolmerSorenson (XX-gate)
                - 'SIS': Square Root I-Swap
                - 'CNOT': CNOT
                - 'ISP': Parameterized I-Swap
                - 'SYC': Sycamore
                - 'SYCP': Parameterized Sycamore

    class functions:
        generate_initial_state(numpy array corresponding to state): this
                generates the appropriate quantum state for the circuits
        two_qubit_rotation(): returns the appropriate rotation for the
                    qubit-qubit coupling
        tfield_rotation(): returns the appropriate transverse field rotation
                    that can be applied
        boundary_rotation(): returns the appropriate boundary field rotation
                    that can be applied
        generate_trotter(time_steps, dt, qreg, qcircuit): generates the
                    trotterization for a given choice of boundary conditions
        next_xx_rotations_pbc(dt, qreg, qcircuit): this function is called in
                    generate trotter to nest the xx rotations with periodic
                    boundary conditions
        next_xx_rotations_obc(dt, qreg, qcircuit): this function is called in
                    generate trotter to nest the xx rotations with open
                    boundary conditions
    '''

    def __init__(self, Jcoup, Gamma, nx, ny, obc=True, xxbasis=True,
                 decompbasis='MS'):
        self.nx, self.ny = nx, ny
        self.efield = Jcoup
        self.mfield = Gamma
        self.obc = obc
        self.xxbasis = xxbasis
        self.decompbases = ['MS', 'SIS', 'CNOT', 'CZ', 'ISP', 'SYC', 'SYCP']
        self.decomps = decompbasis


    def generate_initial_state(self, initial_state):
        circuit = cirq.Circuit()
        qubits = cirq.GridQubits.rect(self.nx, self.ny)
        print('incomplete need to initialize state')
        return


    def generate_initial_state_correlator(self, source, ancilla):
        '''
        a generator to take |psi>|0> -> |psi>|0> + plaquette_op|psi>|1>
        args:
            source (Qubit): the qubit we want to apply the plaquette
                operator on
            ancilla (Qubit): the ancilla to differentiate the states
        '''
        yield cirq.decompose(cirq.H(ancilla))
        # I am assuming a CZ gate here
        yield cirq.decompose(cirq.H(source))
        yield cirq.CZ(ancilla, source)
        yield cirq.decompose(cirq.H(source))



    def two_qubit_rotation(self):
        '''
        this returns the appropriate rotation function from decompositions
        governed by the decomps flag
        '''
        # Assume no decomposition of exp(i s ZZ) gates
        if self.decomps == 'MS':
            if self.xxbasis:
                return XX_native
            else:
                return ZZ_native
        # if we are using the square root iswap CNOT implementation
        if self.decomps == 'SIS':
            if self.xxbasis:
                return decompositions.XX_using_sqrtISWAP
            else:
                return decompositions.ZZ_using_sqrtISWAP
        # if we are using the parameterized iswap rotations
        if self.decomps == 'ISP':
            if self.xxbasis:
                return decompositions.XX_using_param_ISWAP
            else:
                return decompositions.ZZ_using_param_ISWAP
        # if we are using the sycamore CNOT implementation
        if self.decomps == 'SIS':
            if self.xxbasis:
                return decompositions.XX_using_sycamore
            else:
                return decompositions.ZZ_using_sycamore
        # if we are using the parameterized sycamore rotations
        if self.decomps == 'ISP':
            if self.xxbasis:
                return decompositions.XX_using_param_sycamore
            else:
                return decompositions.ZZ_using_param_sycamore
        # if we are using the native CNOT implementation
        if self.decomps == 'CNOT':
            if self.xxbasis:
                return decompositions.XX_using_CNOT
            else:
                return decompositions.ZZ_using_CNOT
        # if we are using the native CNOT implementation
        if self.decomps == 'CZ':
            if self.xxbasis:
                return decompositions.XX_using_CZ
            else:
                return decompositions.ZZ_using_CZ


    def tfield_rotation(self):
        '''
        returns the transverse field rotation
        '''
        if self.xxbasis:
            return cirq.rz
        else:
            return cirq.rx

    def boundary_rotation(self):
        '''
        returns the rotation that is applied to a qubit on the
        boundary
        '''
        if self.xxbasis:
            return cirq.rx
        else:
            return cirq.rz


    def generate_trotter(self, time_steps, dt, qubits):
        '''
        this function produces agenerates a Trotter circuit for the given
        number of time steps
        parameters:
            time_steps (int): the number of time steps
            dt (float): the trotter step size
            qubits (list): the qubits in the circuit
        '''
        circuit = cirq.Circuit()
        # generate the requisit rotation functions
        two_qubit_rot = self.two_qubit_rotation()
        tfield_rot = self.tfield_rotation()
        boundary_rot = self.boundary_rotation()
        # check if we are using periodic boundary conditions
        dtg = -dt * self.mfield
        if not self.obc:
            # iterate over the number of time steps
            for t in range(time_steps):
                # iterate over the sites
                for x in range(self.nx):
                    for y in range(self.ny):
                        circuit += tfield_rot(dtg * 2).on(qubits[x + self.nx * y])
                # apply the electric field operations
                circuit += self.nest_xx_rotations_pbc(dt, qubits,
                                                      two_qubit_rot)
        else:
            for t in range(time_steps):
                # iterate over the sites
                for x in range(self.nx):
                    for y in range(self.ny):
                        scale = 0
                        if (x == 0 or x == self.nx - 1):
                            scale += 1
                        if (y == 0 or y == self.ny - 1):
                            scale += 1
                        circuit += tfield_rot(dtg * 2).on(qubits[x + self.nx * y])
                        # boundary condition effect
                        if scale != 0:
                            circuit += boundary_rot(dt * self.efield * scale * 2).on(qubits[x + self.nx * y])
                # apply the electric field operations
                circuit += self.nest_xx_rotations_obc(dt, qubits,
                                                      two_qubit_rot)
        return circuit


    def nest_xx_rotations_obc(self, dt, qubits, two_q_rot):
        '''
        this function returns a generator that does the efficient nesting for the
        two qubit rotations using periodic boundary conditions
        parameters:
            dt (float): the time step size
            qubits (list): qubit register
            two_q_rot: the two qubit rotation
        '''
        # check if x dimension is even
        for x in range(0, self.nx, 2):
            if x + 1 == self.nx:
                continue
            for y in range(self.ny):
                q1, q2 = qubits[x + self.nx * y], qubits[(x + 1) % self.nx + self.nx * y]
                yield two_q_rot(q1, q2, -dt * self.efield)
        if self.nx % 2 == 0:
            for x in range(1, self.nx - 1, 2):
                if x + 1 == self.nx:
                    continue
                for y in range(self.ny):
                    q1, q2 = qubits[x + self.nx * y], qubits[(x + 1) % self.nx + self.nx * y]
                    yield two_q_rot(q1, q2, -dt * self.efield)
        else:
            for x in range(1, self.nx, 2):
                if x + 1 == self.nx:
                    continue
                for y in range(self.ny):
                    q1, q2 = qubits[x + self.nx * y], qubits[(x + 1) % self.nx + self.nx * y]
                    yield two_q_rot(q1, q2, -dt * self.efield)
        # apply the xx rotations along the other axis
        for x in range(self.nx):
            for y in range(1, self.ny, 2):
                if y + 1 == self.ny:
                    continue
                q1, q2 = qubits[x + self.nx * y], qubits[x + self.nx * ((1 + y) % self.ny)]
                yield two_q_rot(q1, q2, -dt * self.efield)
        if self.ny % 2 == 0:
            for x in range(self.nx):
                for y in range(0, self.ny - 1, 2):
                    if y + 1 == self.ny:
                        continue
                    q1, q2 = qubits[x + self.nx * y], qubits[x + self.nx * ((1 + y) % self.ny)]
                    yield two_q_rot(q1, q2, -dt * self.efield)
        else:
            for x in range(self.nx):
                for y in range(0, self.ny, 2):
                    if y + 1 == self.ny:
                        continue
                    q1, q2 = qubits[x + self.nx * y], qubits[x + self.nx * ((1 + y) % self.ny)]
                    yield two_q_rot(q1, q2, -dt * self.efield)


    def nest_xx_rotations_pbc(self, dt, qubits, two_q_rot):
        '''
        this function does the efficient nesting for the xx rotations
        using periodic boundary conditions
        parameters:
            dt (float): the time step size
            qreg (quantum register): qubit register
            two_q_rot (function): the two qubit rotation
        '''
        # check if x dimension is even
        for x in range(0, self.nx, 2):
            for y in range(self.ny):
                q1, q2 = qubits[x + self.nx * y], qubits[(x + 1) % self.nx + self.nx * y]
                yield two_q_rot(q1, q2, -dt * self.efield)
        if self.nx % 2 == 0:
            for x in range(1, self.nx + 1, 2):
                for y in range(self.ny):
                    q1, q2 = qubits[x + self.nx * y], qubits[(x + 1) % self.nx + self.nx * y]
                    yield two_q_rot(q1, q2, -dt * self.efield)
        else:
            for x in range(1, self.nx, 2):
                for y in range(self.ny):
                    q1, q2 = qubits[x + self.nx * y], qubits[(x + 1) % self.nx + self.nx * y]
                    yield two_q_rot(q1, q2, -dt * self.efield)
        # apply the xx rotations along the other axis
        for x in range(self.nx):
            for y in range(0, self.ny, 2):
                q1, q2 = qubits[x + self.nx * y], qubits[x + self.nx * ((1 + y) % self.ny)]
                yield two_q_rot(q1, q2, -dt * self.efield)
        if self.ny % 2 == 0:
            for x in range(self.nx):
                for y in range(1, self.ny + 1, 2):
                    q1, q2 = qubits[x + self.nx * y], qubits[x + self.nx * ((1 + y) % self.ny)]
                    yield two_q_rot(q1, q2, -dt * self.efield)
        else:
            for x in range(self.nx):
                for y in range(1, self.ny, 2):
                    q1, q2 = qubits[x + self.nx * y], qubits[x + self.nx * ((1 + y) % self.ny)]
                    yield two_q_rot(q1, q2, -dt * self.efield)


class Z2GaugeCirq(object):
    '''
    This object constructs the quantum circuit in the qiskit language for the gauge representation



    class parameters:
        nx, ny: the x and y dimension of the state
        mfield: the spacial plaquette strength (magnetic field)
        efield: the temporal plaquette strength (electric field)
        num_qubits: the number of qubits (links)
        num_plaqs: the number of plaquettes
        obc: a boolean flag indicating whether or not to use open boundary conditions
            default is true
        plaquettes: a list containing the qubits that are associated with a given plaquette
    class functions:
        generate_initial_state(numpy array corresponding to state): this generates the appropriate quantum state
            for the circuits
        generate_trotter(time_steps, dt, qreg, qcircuit): generates the trotterization for a
            given choice of boundary conditions
        oracle_gauss_law(): implement the gauss law constraint to force us onto the physical hilbert space
    '''

    def __init__(self, efield, mfield, nx, ny, obc=True):
        '''
        for now only use OBC (i haven't worked out the pbc implementation)
        '''
        self.efield, self.mfield = efield, mfield
        self.nx, self.ny = nx, ny
        self.obc = obc
        if self.obc:
            self.num_qubits = self.nx * (self.ny - 1) + self.ny * (self.nx - 1)
            self.num_plaqs = (self.nx - 1) * (self.ny - 1)
        # generating an array to indicate which qubits correspond to which links
        array = np.zeros((self.nx, self.ny, 2), dtype='int64') - 1
        counter = 0
        for x in range(self.nx):
            for y in range(self.ny):
                if x < self.nx - 1:
                    array[x, y, 0] = counter
                    counter += 1
                if y < self.ny - 1:
                    array[x, y, 1] = counter
                    counter += 1
        # plaquettes is a list containing tuples which correspond to the qubits
        # which compose the given plaquettes
        self.plaquettes = []
        for x in range(self.nx - 1):
            for y in range(self.ny - 1):
                qubits = (array[x, y, 0], array[x + 1, y, 1],
                          array[x, y + 1, 0], array[x, y, 0])
                self.plaquettes.append(qubits)


    def generate_initial_state(self, initial_state):
        circuit = cirq.Circuit()
        qubits = cirq.GridQubits.rect(self.nx, self.ny)
        print('incomplete need to initialize state')
        return


    def generate_trotter(self, time_steps, dt, qreg):
        '''
        this function generates a Trotter circuit for the the given number of time
        steps
        parameters:
            time_steps (int): the number of time steps
            dt (float): the trotter step size
            qreg (QuantumRegister): the quantum register for the circuit
            qcircuit (QuantumCircuit): the quantum circuit we are operating on
        '''
        circuit = cirq.Circuit()
        for qubit in qreg:
            circuit += cirq.X(qubit) ** (-dt * self.efield)
        ####################################################
        # the plaquette nesting needs to be optimized ######
        ####################################################
        for plaquette in self.plaquettes:
            circuit += cirq.CX(qreg[plaquette[0]], qreg[plaquette[1]])
            circuit += cirq.CX(qreg[plaquette[1]], qreg[plaquette[2]])
            circuit += cirq.CX(qreg[plaquette[2]], qreg[plaquette[3]])
            circuit += cirq.Z(qreg[plaquette[3]]) ** (-dt * self.mfield)
            circuit += cirq.CX(qreg[plaquette[2]], qreg[plaquette[3]])
            circuit += cirq.CX(qreg[plaquette[1]], qreg[plaquette[2]])
            circuit += cirq.CX(qreg[plaquette[0]], qreg[plaquette[1]])
        return circuit

    def oracle_gauss_law(self, qreg):
        print('needs to be implemented')
