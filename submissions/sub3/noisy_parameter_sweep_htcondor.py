import argparse
import os
import numpy as np
import time

from z2_sim.src.QuantumCircuits.Cirq_Code import production
from z2_sim.src.QuantumCircuits.Cirq_Code import io
from z2_sim.src.QuantumCircuits.Cirq_Code.noise.zz_crosstalk_v1 import ZZCrossTalkV1
from z2_sim.src.QuantumCircuits.Cirq_Code.noise.two_local_depol import TwoLocalDepol

import qsimcirq
# assert qsimcirq.qsim_gpu is not None

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('-proc', metavar='proc', type=int, nargs=1,
                    help='htcondor process number')
parser.add_argument('-n', metavar='n', type=int, nargs=1,
                    help='length of spacial dimension')
parser.add_argument('-j', metavar='j', type=float, nargs='+',
                    help='Electric Field Coupling (1/beta)')
parser.add_argument('-dt', metavar='dt', type=float, nargs=1,
                    help='Trotter time step')
parser.add_argument('-tstart', metavar='tstart', type=int, nargs=1,
                    help='Trotter step startpoint')
parser.add_argument('-tstop', metavar='tstop', type=int, nargs=1,
                    help='Trotter step stop-point')
parser.add_argument('-zeta', metavar='zeta', type=int, nargs='+',
                    help='Unitary crosstalk noise parameter')
parser.add_argument('-eps', metavar='eps', type=float, nargs='+',
                    help='Depolarizing noise parameter')
parser.add_argument('-r', metavar='r', type=int, nargs=1,
                    help='Number of trajectories to simulate.')
parser.add_argument('-dest', metavar='dest', type=str, nargs=1,
                    help='Directory to save results.')
args = parser.parse_args()

# Subroutine
proc = args.proc[0]
n = args.n[0]
j_sweep = args.j
dt = args.dt[0]
tstart = args.tstart[0]
tstop = args.tstop[0]
zeta_sweep = args.zeta
eps_sweep = args.eps
n_trajectories = args.r[0]
dest = args.dest[0]

# Pre-compose noise models
target_gate="SIS"
GATE_DURATION = 1e-8
zeta_model_sweep = [
    ZZCrossTalkV1(zeta, target_gate, gate_duration=GATE_DURATION, sampled_error=False) for zeta in zeta_sweep
]

eps_model_sweep = [
    TwoLocalDepol(err_1q=eps / 10, err_2q=eps, sampled_error=False, verbose=False) for eps in eps_sweep
]

# Initialize simulators and noise model
N_FUSE = 4
qsim_options = qsimcirq.QSimOptions(
    max_fused_gate_size=N_FUSE,
    ev_noisy_repetitions=n_trajectories,
    use_gpu=True,
    gpu_sim_threads=256,
    gpu_state_threads=512,
    gpu_data_blocks=16,
    verbosity=0,
    denormals_are_zeros=True,
)

for jcoup in j_sweep:
    for i, zeta in enumerate(zeta_sweep):
        for j, eps in enumerate(eps_sweep):
            zeta_model = zeta_model_sweep[i]
            eps_noise_model = eps_model_sweep[j]
            t0 = time.time()
            # TODO: if eps == 0 we need to manually implement noiseless simulation
            if eps < 1e-9:
                # Reroute to a unitary simulation with intermediate state vectors.
                out = production.compute_obs_with_intermediate_state_vector(
                    n=n,
                    trotter_steps=tstop - tstart,
                    jcoup=jcoup,
                    dt=dt,
                    all_observables=True,
                    qsim_options=dict(t=8, f=4, g=False), # GPU not necessary
                    decompbasis=target_gate,
                    obc=True,
                    noise_models=[zeta_model],
                )
            else:
                out = production.compute_noisy_obs(
                    n=n,
                    trotter_start=tstart,
                    trotter_stop=tstop,
                    jcoup=jcoup,
                    dt=dt,
                    all_observables=True,
                    qsim_options=qsim_options,
                    noise_models=[zeta_model, eps_noise_model],
                    obc=True,
                    decompbasis=target_gate,
                )
            delta_time = time.time() - t0
            print(f"size {n}, one pass, trotter interval=({tstart}, {tstop}), {n_trajectories} trajectories")
            print(f"\t{delta_time}")
            # This is to expedite file transfers from the submission nodes, since the .sub
            # doesn't know the output file name
            fout = io.make_noisy_htcondor_run_fname(proc, n, jcoup, dt, tstart, tstop, zeta, eps, n_trajectories)
            np.save(os.path.join(dest, fout), out)