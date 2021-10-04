import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from multiprocessing import Pool


class Z2GaugeNumpy():
    '''
    class to manage constructing states and running states and for z_2 dual ising gauge using
    numpy simulations

    parameters:
        - nx (int): the x dimension of the lattice
        - ny (int): the y dimension of the lattice
        - nt (int): the temporal extent of the monte carlo lattice we are looking at
        - dim (int): the size of the Hilbert space
        - jcoup (float): the electric field strength
        - gamma (float): the magnetic field strength
        - dt (float): the Trotter step size
        - obc (boolean): whether or not open boundary conditions (obc) are being used
                        by default this is true
        - xxbasis (boolean): a flag on whether to construct the electric field operators
                            being off diagonal (xx) or diagonal (zz)


    functions
        - construct_operators(): constructs the foundational operators necessary for the
                                    the trotter operators and hamiltonian
        - construct_hamiltonian(): constructs the hamiltonian for the z2 gauge theory
        - constructtrotteroperators(): constructs the trotter operators
        - construct_exciteops(): constructs the elementary site operators to create a particle
                                 excitation (in the deconfined basis)
        - load_monte_carlo(num): loads the num^th monte carlo configuration
        - construct_excitation_single_mode(num, kmom): construct the even / odd parity states for single excitation on configuration num, with momentum kmom
        - construct_excitation_two_particle(num, kmom, pmom): construct the even / odd parity states for two excitation on configuration num, with momentum kmom and pmom
        - time_evolve_monte_carlo_single_mode(args): run a single mode monte carlo run
    '''

    def __init__(self, nx, ny, jcoup, gamma, dt, obc=True, xxbasis=True):
        self.nx = nx
        self.ny = ny
        self.nt = 96
        self.dim = 2 ** (nx * ny)
        self.jcoup = jcoup
        self.gamma = gamma
        self.obc = obc
        self.dt = dt
        self.xxbasis = xxbasis
        self.construct_operators()
        self.constructtrotteroperators()
#         self.construct_exciteops()
#         self.construct_hamiltonian()


    def construct_hamiltonian(self):
        '''
        construct the hamiltonian and extract the ground state energy
        '''
        ham = 0
        if self.xxbasis:
            # iterate through the
            for op in self.zoperators:
                ham -= self.gamma * op
            for op in self.xxoperators:
                ham -= self.jcoup * op
            if self.obc:
                for i in range(self.nx):
                    for j in range(self.ny):
                        scale = 0
                        if i == 0 or i == self.nx - 1:
                            scale += 1
                        if j == 0 or j == self.ny - 1:
                            scale += 1
                        ham -= scale * self.jcoup * self.xoperators[i + self.nx * j]
        else:
            # iterate through the
            for op in self.xoperators:
                ham -= self.gamma * op
            for op in self.zzoperators:
                ham -= self.jcoup * op
            if self.obc:
                for i in range(self.nx):
                    for j in range(self.ny):
                        scale = 0
                        if i == 0 or i == self.nx - 1:
                            scale += 1
                        if j == 0 or j == self.ny - 1:
                            scale += 1
                        ham -= scale * self.jcoup * self.zoperators[i + self.nx * j]
        self.ham = ham


    def construct_operators(self):
        '''
        construct the x, xx, and z operators
        '''
        self.zoperators = []
        self.xoperators = []
        self.xxoperators = []
        self.zzoperators = []
        xop = sp.bsr_matrix([[0, 1], [1, 0]], dtype='complex128')
        zop = sp.bsr_matrix([[1, 0], [0, -1]], dtype='complex128')
        self.magop = 0
        nx, ny = self.nx, self.ny
        for i in range(self.nx * self.ny):
            id1 = sp.identity(2**i, dtype='complex128')
            id2 = sp.identity(2**(self.nx * self.ny - 1 - i), dtype='complex128')
            x_i = sp.kron(id1, sp.kron(xop, id2))
            z_i = sp.kron(id1, sp.kron(zop, id2))
            # print(z_i.shape)
            self.zoperators.append(z_i)
            self.magop += z_i / 8
            self.xoperators.append(x_i)


        for x in range(self.nx):
            for y in range(self.ny):
                index1 = x + nx * y
                index2 = (x + 1) % nx + y * nx
                if not (self.obc and x + 1 == nx):
                    # print(index1, index2)
                    if self.xxbasis:
                        xxop = self.xoperators[index1].dot(self.xoperators[index2])
                        self.xxoperators.append(xxop)
                    else:
                        zzop = self.zoperators[index1].dot(self.zoperators[index2])
                        self.zzoperators.append(zzop)
                index3 = x + ((y + 1) % nx) * ny
                if not (self.obc and y + 1 == ny):
                    # print(index1, index3)
                    if self.xxbasis:
                        xxop = self.xoperators[index1].dot(self.xoperators[index3])
                        self.xxoperators.append(xxop)
                    else:
                        zzop = self.zoperators[index1].dot(self.zoperators[index3])
                        self.zzoperators.append(zzop)


    def constructtrotteroperators(self):
        '''
        construct the trotterizations operators
        '''
        self.plaqops = []
        self.boundops = []
        self.linkops = []
        identity = sp.identity(2**(self.nx * self.ny), dtype='complex128')
        # electric field is off diagonal
        if self.xxbasis:
            # magnetic field rotations
            for zop in self.zoperators:
                op = identity * np.cos(self.gamma * self.dt) - zop * np.sin(self.gamma * self.dt) * 1.0j
                self.plaqops.append(op)
            # electric field rotations
            for xxop in self.xxoperators:
                op = identity * np.cos(self.dt * self.jcoup) - xxop * np.sin(self.dt * self.jcoup) * 1.0j
                self.linkops.append(op)
            if self.obc:
                for i in range(self.nx * self.ny):
                    x = i % self.nx
                    y = (i // self.nx) % self.ny
                    scale = 0
                    if x == 0 or x == self.nx - 1:
                        scale += 1
                    if y == 0 or y == self.ny - 1:
                        scale += 1
                    if scale > 0:
                        cos = np.cos(self.dt * self.jcoup * scale)
                        sin = 1.0j * np.sin(self.dt * self.jcoup * scale)
                        self.boundops.append(identity * cos - self.xoperators[i] * sin)
        # electric field is diagonal
        else:
            # magnetic field operators
            for xop in self.xoperators:
                op = identity * np.cos(self.gamma * self.dt) - xop * np.sin(self.gamma * self.dt) * 1.0j
                self.plaqops.append(op)
            # electric field operators
            for zzop in self.zzoperators:
                op = identity * np.cos(self.dt * self.jcoup) - zzop * np.sin(self.dt * self.jcoup) * 1.0j
                self.linkops.append(op)
            # boundary electric field operators
            if self.obc:
                for i in range(self.nx * self.ny):
                    x = i % self.nx
                    y = (i // self.nx) % self.ny
                    scale = 0
                    if x == 0 or x == self.nx - 1:
                        scale += 1
                    if y == 0 or y == self.ny - 1:
                        scale += 1
                    if scale > 0:
                        cos = np.cos(self.dt * self.jcoup * scale)
                        sin = 1.0j * np.sin(self.dt * self.jcoup * scale)
                        self.boundops.append(identity * cos - self.zoperators[i] * sin)


    def measure_wilson_loops(self, timesteps, naivegs=True):
        '''
        perform a Trotter evolution to measure the correlator of two spacial wilson loops

        Args:
            timesteps (int): number of time steps we want to take
            naivegs (bool): whether to use a naive groundstate or not
        '''
        # use the naive ground state
        if not naivegs:
            new1 = sp.linalg.eigsh(self.ham, which='SA', k=1)[1]
        else:
            new1 = np.reshape(np.array([1 if i == 0 else 0 for i in range(2 ** (self.nx * self.ny))],
                                       dtype='complex128'),
                              (1, 2 ** (self.nx * self.ny))).transpose()
        if self.xxbasis:
            new2 = new1#self.zoperators[0].dot(new1)
        else:
            new2 = self.xoperators[0].dot(new1)
        # object to contain the correlator values
        correlator = np.zeros((self.nx, self.ny, timesteps), dtype='complex128')
        for i in range(timesteps):
            for op in self.plaqops:
                new1 = op.dot(new1)
                new2 = op.dot(new2)
            for op in self.boundops:
                new1 = op.dot(new1)
                new2 = op.dot(new2)
            for op in self.linkops:
                new1 = op.dot(new1)
                new2 = op.dot(new2)
            for j in range(self.nx):
                for k in range(self.ny):
                    if self.xxbasis:
                        vec2 = self.zoperators[j * self.ny + k].dot(new1)
                    else:
                        vec2 = self.xoperators[j * self.ny + k].dot(new1)
                    val = vec2.conjugate().transpose().dot(new2)
                    correlator[j, k, i] += val[0, 0]
        return correlator


    def construct_exciteops(self):
        '''
        construct the excitation operators
        '''
        self.excites = {}
        excite_op = sp.bsr_matrix([[0, 0], [1, 0]], dtype='complex128')
        for i in range(self.nx * self.ny):
            operator = sp.kron(sp.identity(2 ** i, dtype='complex128'),
                               sp.kron(excite_op,
                                       sp.identity(2 ** (self.nx * self.ny - i - 1),
                                                   dtype='complex128')))
            self.excites[(i // self.nx, i % self.nx)] = operator



    def load_monte_carlo(self, num):
        '''
        loads the monte carlo configuration and returns the euclidean lattice
        '''
        nt = self.nt
        nx = self.nx
        ny = self.ny
        filename = '../configurations/configuration{}'.format(num)
        filename += 'j=0_3ht=1_0'.format(nx)
        filename += 'ns_{}_nt_{}'.format(nx, nt)
        filename += '.csv'
        config = np.genfromtxt(filename, delimiter=',')
        lattice = np.zeros((nx, ny, nt))
        for t in range(nt):
            for i in range(nx * ny):
                x = i // nx
                y = i % nx
                lattice[x, y, t] = config[t, i]
        return lattice


    def construct_excitation_single_mode(self, num, kmom):
        '''
        construct a single momentum mode excitation on the lattice
        '''
        nx, ny = self.nx, self.ny
        kx, ky = kmom[0], kmom[1]
        dual = self.convert_to_dual(num)
        dim = self.dim
        state1 = sp.dok_matrix((1, dim), dtype='complex128')
        state2 = sp.dok_matrix((1, dim), dtype='complex128')
        sym = 0
        asym = 0
        lat1 = np.array(dual[:, :, 3], copy=True)
        lat2 = np.array(dual[:, :, -3], copy=True)
        # initialize the state
        useflag = True
        index1, index2 = 0, 0 #2 ** (nx * ny), 2 ** (nx * ny)
        for i in range(nx * ny):
            if lat1[i // nx, i % nx] == -1:
                index1 += 2 ** i
            if lat2[i // nx, i % nx] == 1:
                index2 += 2 ** i
        state1[0, index1] = 1.0
        state2[0, index2] = 1.0

        # convert to a compressed sparse column matrix
        state1 = state1.tocsc().transpose()
        state2 = state2.tocsc().transpose()
        # generate the excitation operator
        excite_op = sp.bsr_matrix((2 ** (nx * ny), 2 ** (nx * ny)),
                                  dtype='complex128')
        for x in range(nx):
            for y in range(ny):
                excite_op += self.excites[(x, y)] * np.exp(-2.0j * np.pi * (x * kx / nx + y * ky / ny))
        # apply the excitation operator
        state1 = excite_op.dot(state1)
        state2 = excite_op.dot(state2)
        counter1, counter2 = 0, 0
        # find the norms of the states
        norm1 = np.sqrt((state1.conjugate().transpose().dot(state1)).todense()[0, 0])
        norm2 = np.sqrt((state2.conjugate().transpose().dot(state2)).todense()[0, 0])
        if norm1 != 0:
            counter1 += 1
            state1 /= norm1
        if norm2 != 0:
            counter2 += 1
            state2 /= norm2
        sym = state1 + state2
        asym
        if counter1 != 0 and counter2 != 0:
            # make symmetric and anti-symmetric states
            sym = state1 + state2
            asym = state1 - state2
            # construct asymflag
            norm = np.sqrt((sym.conjugate().transpose().dot(sym)).todok()[0, 0])
            sym /= norm
            isasym = True
            # check if state1 and state2 are the same
            norm2 = np.sqrt((asym.dot(asym.conjugate().transpose())).todok()[0, 0])
            if  norm2 == 0:
                isasym = False
                return sym, np.array([[0]])
            else:
                asym /= np.sqrt(norm2)
                return sym, asym
        else:
            return np.array([[0]]), np.array([[0]])

    def construct_excitation_two_particle(self, num, kmom, lmom):
        '''
        constructs the two particle plane waves
        '''
        nx, ny = self.nx, self.ny
        kx, ky = kmom[0], kmom[1]
        qx, qy = lmom[0], lmom[1]
        dual = self.convert_to_dual(num)
        dim = self.dim
        state1 = sp.dok_matrix((1, dim), dtype='complex128')
        state2 = sp.dok_matrix((1, dim), dtype='complex128')
        sym = 0
        asym = 0

        lat1 = np.array(dual[:, :, 4], copy=True)
        lat2 = np.array(dual[:, :, 43], copy=True)
        useflag = True
        # prepare the momentum state
        counter1, counter2 = 0, 0
        for vx in range(nx):
            for vy in range(ny):
                for wx in range(nx):
                    for wy in range(ny):
                        # don't do anything if both excitations are at
                        # the same spot
                        if vx == wx and vy == wy:
                            continue
                        else:
                            copy1 = np.array(lat1)
                            copy2 = np.array(lat2)
                            # check if we can make an excitation
                            if lat1[vx, vy] == 1 and lat1[wx, wy] == 1:
                                copy1[vx, vy] = -1
                                copy1[wx, wy] = -1
                                counter1 += 1
                            if lat2[vx, vy] == 1 and lat2[wx, wy] == 1:
                                copy2[vx, vy] = -1
                                copy2[wx, wy] = -1
                                counter2 += 1
                            index1 = 0
                            index2 = 0
                            for i in range(nx * ny):
                                if copy1[i // nx, i % nx] == -1:
                                    index1 += 2 ** i
                                if copy2[i // nx, i % nx] == -1:
                                    index2 += 2 ** i
                            t1 = np.exp(-2.0j * np.pi * vx * kx / nx + 1.0j * np.pi / nx)
                            t2 = np.exp(-2.0j * np.pi * vy * ky / ny + 1.0j * np.pi / ny)
                            t3 = np.exp(-2.0j * np.pi * wx * qx / nx + 1.0j * np.pi / nx)
                            t4 = np.exp(-2.0j * np.pi * wy * qy / ny + 1.0j * np.pi / ny)
                            t5 = np.exp(-2.0j * np.pi * wx * kx / nx + 1.0j * np.pi / nx)
                            t6 = np.exp(-2.0j * np.pi * wy * ky / ny + 1.0j * np.pi / ny)
                            t7 = np.exp(-2.0j * np.pi * vx * qx / nx + 1.0j * np.pi / nx)
                            t8 = np.exp(-2.0j * np.pi * vy * qy / ny + 1.0j * np.pi / ny)
                            term = t1 * t2 * t3 * t4 - t5 * t6 * t7 * t8
                            state1[0, index1] += term
                            state2[0, index2] += term
        state1 = state1.tocsc()
        state2 = state2.tocsc()
        # print(state1)
        if counter1 != 0 and counter2 != 0:
            # make symmetric and anti-symmetric states
            sym = state1 + state2
            asym = state1 - state2
            # construct asymflag
            norm = np.sqrt(sym.dot(sym.conjugate().transpose())[0, 0])
            #print(norm)
            sym /= norm
            isasym = True
            # print(sym.dot(sym.conjugate().transpose()))
            # check if state1 and state2 are the same
            norm2 = np.sqrt(asym.dot(asym.conjugate().transpose())[0, 0])
            if  norm2 == 0:
                isasym = False
                return sym, np.array([[0]])
            else:
                asym /= np.sqrt(norm2)
                return sym, asym
        else:
            return np.array([[0]]), np.array([[0]])

    def time_evolve_monte_carlo_single_mode(self, args):
        '''
        carry out the time evolution of a monte carlo state this is designed to be run in parallel
        args is a tuple of the form:
        (time_steps (int), config_number (int), momentum_state (int, int))
        '''
        tsteps = args[0]
        num = args[1]
        k = args[2]
        magsmontecarlo = np.zeros(tsteps, dtype='complex128')
        symevo, asymevo = self.construct_excitation_single_mode(num, k)
        isymevo = sp.csr_matrix(symevo, copy=True).conjugate().transpose()
        iasymevo = sp.csr_matrix(asymevo, copy=True).conjugate().transpose()
        print('configuration {}'.format(num))
        if symevo.shape[0] == 1:
            return (magsmontecarlo, 0)
        else:
            for i in range(tsteps):
                val = 0
                for op in self.plaqops:
                    symevo = op.dot(symevo)
                for op in self.linkops:
                    symevo = op.dot(symevo)
                if asymevo.shape[0] != 1:
                    for op in self.plaqops:
                        asymevo = op.dot(asymevo)
                    for op in self.linkops:
                        asymevo = op.dot(asymevo)
                    val += (isymevo.dot(symevo))[0, 0]
                    val -= (iasymevo.dot(asymevo))[0, 0]
                    magsmontecarlo[i] += val / 2
                else:
                    val += (isymevo.dot(symevo))[0, 0]
                    magsmontecarlo[i] += val
            return (magsmontecarlo, 1)


