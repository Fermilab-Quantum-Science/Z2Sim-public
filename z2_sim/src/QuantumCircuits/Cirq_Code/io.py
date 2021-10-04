"""io.py - input/output management for saving and loading data."""
import os
from typing import Optional, Sequence, Tuple

import numpy as np

def make_sweep_fname(tag: str, n: int, obc: bool, noisy: Optional[str] = False):
    """Return a filename for noiseless parameter sweeps outputs.

    For cross-compatibility I'm opting to store sweeps as two files: The first
    file contains the array of sweeped parameters
    Args:
        tag: Some identifier string for this simulation result.
        n: Grid size simulated.
        obc: Boolean indicating Open Boundary Conditions (OBC) if True,
            otherwise Periodic Boundary Conditions (PBC).
    """

    bcstr = "obc" if obc else "pbc"
    noisestr = "noisy" if noisy else "noiseless"
    return f"{noisestr}_parameter_sweep_{tag}_results_n{n}_{bcstr}.npy"


def make_physical_params_sweep_fname(tag: str, param_str: str):
    """Return a filename for noiseless parameter sweeps inputs.

    Args:
        tag: Some identifier string for this simulation result.
    """
    return f"physical_parameter_sweep_{tag}_{param_str}.npy"


def make_noisy_params_sweep_fname(tag: str, param_str: str):
    """Return a filename for noiseless parameter sweeps inputs.

    Args:
        tag: Some identifier string for this simulation result.
    """
    return f"noisy_parameter_sweep_{tag}_{param_str}.npy"


def make_noisy_htcondor_run_fname(
    proc: int,
    n: int,
    j: float,
    dt: float,
    tstart: int,
    tstop: int,
    zeta: int,
    eps: float,
    r: int,
    ):
    """Return a filename for noiseless parameter sweeps inputs.

    Args:
        TODO
    """

    jrd = np.round(j, decimals=6)
    dtrd = np.round(dt, decimals=6)
    epsrd = np.round(eps, decimals=6)

    # Legacy:
    return (f"results_proc{proc}_n{n}_j{jrd}_dt{dtrd}_tstart{tstart}"
                + f"_tstop{tstop}_zeta{zeta}_eps{epsrd}.npy")
    # return (f"results_proc{proc}_n{n}_j{jrd}_dt{dtrd}_tstart{tstart}"
    #             + f"_tstop{tstop}_zeta{zeta}_eps{epsrd}_r{r}.npy")


class CondorCollector:
    """Manage the ouput files of a set of a multinode htcondor run."""
    def __init__(
        self,
        path: str,
        nprocs: int,
        n: int,
        j_sweep: Sequence[float],
        dt_sweep: Sequence[float],
        trotter_intervals: Sequence[Tuple[int]],
        zeta_sweep: Sequence[int],
        eps_sweep: Sequence[float],
        r: int,
    ):
        """
        Args:

        """

        self.path = path
        self.nprocs = nprocs
        self.n = n
        self.j_sweep = j_sweep
        self.dt_sweep = dt_sweep
        self.trotter_intervals = trotter_intervals
        self.zeta_sweep = zeta_sweep
        self.eps_sweep = eps_sweep
        self.r = r

    def collect_htcondor_outputs(self):
        """Collect all htcondor outputs into a single .npy.

        Returns:
            An array of the shape

                `(nprocs, len(j_sweep), len(dt_sweep), len(zeta_sweep), len(eps_sweep), trotter_steps, n, n)`

        compute_noisy_obs
                (trotter_stop - trotter_start, n, n)
        """

        # The trotter steps of many runs will be "collapsed" into a single run
        n_trotter = self.trotter_intervals[-1][1] - self.trotter_intervals[0][0]
        out = np.zeros((
            self.nprocs,
            len(self.j_sweep),
            len(self.dt_sweep),
            len(self.zeta_sweep),
            len(self.eps_sweep),
            n_trotter,
            self.n,
            self.n,
        ), dtype=complex)
        for p in range(self.nprocs):
            for i, jcoup in enumerate(self.j_sweep):
                for j, dt in enumerate(self.dt_sweep):
                    for k, zeta in enumerate(self.zeta_sweep):
                        for ell, eps in enumerate(self.eps_sweep):
                            # print(i, " ", j, " ", k, " ", ell, " ")
                            # print(ell + 6*k + 36 * i)

                            temp = []
                            for (tstart, tstop) in self.trotter_intervals:
                                target = make_noisy_htcondor_run_fname(
                                    proc=p,
                                    n=self.n,
                                    j=jcoup,
                                    dt=dt,
                                    tstart=tstart,
                                    tstop=tstop,
                                    zeta=zeta,
                                    eps=eps,
                                    r=self.r
                                    )
                                x = np.load(os.path.join(self.path, target))

                                # try:
                                #     x = np.load(os.path.join(self.path, target))
                                # except FileNotFoundError:
                                #     print("whoops! missing jcoup={}, zeta={}, eps={}".format(jcoup, zeta, eps))
                                temp.append(x)
                            out[p, i, j, k, ell] = np.vstack(temp)

        return out