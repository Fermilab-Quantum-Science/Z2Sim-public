import argparse
import os
import numpy as np
import time

from z2_sim.src.QuantumCircuits.Cirq_Code import production
from z2_sim.src.QuantumCircuits.Cirq_Code import io
from z2_sim.src.QuantumCircuits.Cirq_Code.noise.zz_crosstalk_v1 import ZZCrossTalkV1
from z2_sim.src.QuantumCircuits.Cirq_Code.noise.two_local_depol import TwoLocalDepol

import qsimcirq

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('-gpu', metavar='gpu', type=int, nargs=1,
                    help='boolean indicating GPU')
parser.add_argument('-resource', metavar='resource', type=str, nargs=1,
                    help='string tag describing this device')
args = parser.parse_args()

# Fix noise and physical parameters
#
gpu = args.gpu[0]
if gpu:
    assert qsimcirq.qsim_gpu is not None

resource = args.resource[0]
jcoup = 0.8
dt = 0.25
tstart = 1
tstop = 10
zeta = 10000
# max this out since technically thats the worst-case
eps = 0.005

# Pre-compose noise models
target_gate="SIS"
GATE_DURATION = 1e-8
zeta_model = ZZCrossTalkV1(zeta, target_gate, gate_duration=GATE_DURATION, sampled_error=False)

eps_noise_model = TwoLocalDepol(err_1q=eps / 10, err_2q=eps, sampled_error=False, verbose=False)

# qsim_options = qsim_options={'r': trajectory_count, 't': N_THREADS, 'f': N_FUSE}
# Initialize simulators and noise model
N_FUSE = 4
# # # # # # # # # # # # # # # # # # # # # # # # #
# MODIFY THIS BEFORE RUNNING LOCALLY
N_THREADS = 32
# # # # # # # # # # # # # # # # # # # # # # # # #

trajectory_sweep = [500]
trajectory_trials = [1]
grid_sweep = [5]

# First axis is for trajectories
# Second axis is for grid sizes = sqrt(total qubits)
times = np.zeros((len(trajectory_sweep), len(grid_sweep)))

for xx, (n_trajectories, ntrials) in enumerate(zip(trajectory_sweep, trajectory_trials)):
    if gpu == 0:
        qsim_options = qsimcirq.QSimOptions(
            max_fused_gate_size = N_FUSE,
            cpu_threads = N_THREADS,
            ev_noisy_repetitions = n_trajectories,
            denormals_are_zeros = True,
        )
    elif gpu == 1:
        qsim_options = qsimcirq.QSimOptions(
            max_fused_gate_size = N_FUSE,
            use_gpu=True,
            ev_noisy_repetitions = n_trajectories,
            denormals_are_zeros = True,
        )
    for yy, n in enumerate(grid_sweep):
        temp = []
        # Expect much bigger fluctuations for lower rep count.
        for _ in range(ntrials):
            t0 = time.time()
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
            temp.append(delta_time)
            print(f"size {n}, one pass, trotter interval=({tstart}, {tstop}), {n_trajectories} trajectories")
            print(f"\t{delta_time}")

        times[xx, yy] = min(temp)

np.save(f"benchmark_gpu{gpu}_{resource}_tstart{tstart}.npy", times)