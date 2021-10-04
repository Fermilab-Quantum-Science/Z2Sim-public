import os
import shutil

import numpy as np

from z2_sim.src.QuantumCircuits.Cirq_Code import io


def test_condor_collector():
    """Condor collector should appropriately grab and combine outputs."""

    nprocs = 2
    n = 4
    tstart, tstop = (1, 11)
    n_trajectories = 33
    dt = 0.25
    j_sweep = [1.0, 2, 3.22]
    zeta_sweep = [44, 55]
    eps_sweep = [.001, .0005]

    dest = "./temp"
    os.mkdir(dest)
    for proc in range(nprocs):
        for jcoup in j_sweep:
            for zeta in zeta_sweep:
                for eps in eps_sweep:
                    dummy = np.ones((tstop - tstart, n, n))
                    fout = io.make_noisy_htcondor_run_fname(proc, n, jcoup, dt, tstart, tstop, zeta, eps, n_trajectories)
                    np.save(os.path.join(dest, fout), dummy)

    collector = io.CondorCollector(
        path=dest,
        nprocs=nprocs,
        n=n,
        j_sweep=j_sweep,
        dt_sweep=[dt],
        trotter_intervals=[(tstart, tstop)],
        zeta_sweep=zeta_sweep,
        eps_sweep=eps_sweep,
        r=n_trajectories,
    )
    out = collector.collect_htcondor_outputs()
    assert out.shape == (nprocs, len(j_sweep), 1, len(zeta_sweep), len(eps_sweep), tstop - tstart, n, n)
    shutil.rmtree(dest)


if __name__ == "__main__":
    test_condor_collector()