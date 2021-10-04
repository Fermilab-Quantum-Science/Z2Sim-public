# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 11:12:29 2021

@author: gusta
"""


import argparse
# @Erik let me know if there's problems with this style import?
from z2_sim.src.NumericsPython import Z2Gauge
# import Z2Gauge
import numpy as np
import scipy.sparse as sp


def z2_dual_correlater_numeric(nspace, Jcoupling, Gamma, deltat, ntime, obc, localgroundstate, exgs=False):
    lattice = Z2Gauge.Z2GaugeNumpy(nspace, nspace,
                                   Jcoupling, Gamma,
                                   deltat, obc=obc,
                                   xxbasis=False)
    lattice.construct_hamiltonian()
    if exgs:
        gs = sp.csr_matrix(eigsys[1][:, 0]).transpose()
        if nspace % 2 == 1:
            index = (nspace * nspace - 1) // 2 + 1
            e1 = lattice.xoperators[index].dot(gs)
#         print(gs.shape, e1.shape)
    # here we need to generate the initial state
    elif localgroundstate:
        eigsys = np.linalg.eigh(-np.array([[4 * Jcoupling,
                                            Gamma],
                                           [Gamma,
                                            -4 * Jcoupling]]))
        gs_plaq = eigsys[1][:, 0]
        gs = np.kron(np.kron(np.kron(np.kron(gs_plaq, gs_plaq),
                                     np.kron(gs_plaq, gs_plaq)),
                             np.kron(np.kron(gs_plaq, gs_plaq),
                                     np.kron(gs_plaq, gs_plaq))),
                     np.kron(np.kron(np.kron(gs_plaq, gs_plaq),
                                     np.kron(gs_plaq, gs_plaq)),
                             np.kron(np.kron(gs_plaq, gs_plaq),
                                     np.kron(gs_plaq, gs_plaq))))
        gs = sp.bsr_matrix(gs)
        # generate the excited state for a 1 x 1 plaquette
        if nspace % 2 == 1:
            index = (nspace * nspace - 1) // 2 + 1
            e1 = lattice.xoperators[index].dot(gs)
        else:
            index = (nspace ** 2) // 2
            e1 = lattice.xoperators[index].dot(gs)
    else:
        gs = sp.dok_matrix((2 ** (nspace * nspace), 1))
        gs[0, 0] = 1
        gs.tocsc()
        if nspace % 2 == 1:
            index = (nspace * nspace - 1) // 2
            e1 = lattice.xoperators[index].dot(gs)
        else:
            index = (nspace ** 2) // 2 + 1
            e1 = lattice.xoperators[index].dot(gs)
    # here we for the time evolution
    correlator = np.zeros((nspace, nspace,
                           ntime + 1),
                          dtype='complex128')
    # Trotter evolution
    for index in range(nspace ** 2):
        x = index % nspace
        y = index // nspace
        op = lattice.xoperators[index]
        val = e1.conjugate().transpose().dot(op).dot(gs)
#         correlator[x, y, 0] += val
#         print(val)
        correlator[x, y, 0] += val[0,0]
    for t in range(ntime):
        # apply the z-operators
        # apply the x-type operators
        for op in lattice.plaqops:
            e1 = op.dot(e1)
            gs = op.dot(gs)
        for op in lattice.boundops:
            e1 = op.dot(e1)
            gs = op.dot(gs)
        for op in lattice.linkops:
            e1 = op.dot(e1)
            gs = op.dot(gs)
        # now we calculate the correlator values
        for index in range(nspace ** 2):
            x = index % nspace
            y = index // nspace
            op = lattice.xoperators[index]
            val = e1.conjugate().transpose().dot(op).dot(gs)
#             print(val)
            correlator[x, y, t + 1] += val[0,0]
#             correlator[x, y, t + 1] += val

    return correlator


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('nspace', metavar='NS', type=int, nargs=1,
                        help='length of spacial dimension')
    parser.add_argument('ntime', metavar='NT', type=int, nargs=1,
                        help='number of Trotter steps')
    parser.add_argument('Jcoupling', metavar='J', type=float, nargs=1,
                        help='Electric Field Coupling (beta)')
    parser.add_argument('Gamma', metavar='G', type=float, nargs=1,
                        help='Magnetic Field Coupling (1/beta)')
    parser.add_argument('deltat', metavar='dt', type=float, nargs=1,
                        help='Trotter time step')
    parser.add_argument('-obc', '--obc', action='store_true',
                        default=False)
    parser.add_argument('-lgs', '--localgroundstate', action='store_true',
                        default=False)

    args = parser.parse_args()

    # Subroutine
    correlator = z2_dual_correlater_numeric(args.nspace[0], args.Jcoupling[0], args.Gamma[0], args.deltat[0], args.ntime[0], args.periodic, args.localgroundstate)

    # save the correlators as a numpy array
    filename = 'CorrelatorJ={}Gamma={}'.format(args.Jcoupling[0], args.Gamma[0])
    filename += 'dt={}ns={}nt={}'.format(args.deltat[0], args.nspace[0], args.ntime[0])
    filename = filename.replace('.', '_')
    filename += '.npy'
    np.save(filename, correlator)
