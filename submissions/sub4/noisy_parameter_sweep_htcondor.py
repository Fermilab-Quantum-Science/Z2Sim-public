"""noisy_parameter_sweep_htcondor.py

NOTE: Jobs parameters (sweep values, repetitions) must be modified directly in
this script. This script is configured to be dispatched based on a single
int $(Process) passed in by the htcondor submission.

"""
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
parser.add_argument('-dest', metavar='dest', type=str, nargs=1,
                    help='Directory to save results.')
args = parser.parse_args()

dest = args.dest[0]
### DISPATCHER ###
proc = args.proc[0]
##################

n = 4
tstart = 1
tstop = 51
n_trajectories = 1000

# Hardcoded physical parameters based on Hank's input
dt = 0.25

#### Hardcoded TABLE of noise parameter sweep ####
j_sweep = [0.714285, 0.625, .555556]
zeta_sweep = [0, 150000, 300000, 450000, 600000, 750000]
eps_sweep = [0, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003]
##################################################
j_zeta_eps_table = []
for j in j_sweep:
    for zeta in zeta_sweep:
        for eps in eps_sweep:
            j_zeta_eps_table.append((j, zeta, eps))

# Now dispatch.
jcoup, zeta, eps = j_zeta_eps_table[proc]
print(f"PROC {proc}: (zeta, epsilon)=({zeta}, {eps})")
target_gate="SIS"
GATE_DURATION = 1e-8

# Pre-compose noise models
zeta_model = ZZCrossTalkV1(zeta, target_gate, gate_duration=GATE_DURATION, sampled_error=False)
eps_noise_model = TwoLocalDepol(err_1q=eps / 10, err_2q=eps, sampled_error=False, verbose=False)


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


t0 = time.time()
# TODO: if eps == 0 we need to manually implement noiseless simulation
if eps < 1e-9:
    print("Dispatching to noiseless simulator for eps=0")
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
    print("Dispatching to triangle simulator for eps!=0")
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

# !!!!!!!!
# Proc is hardcoded to `0` in the file output for forwards compatibility
# with the CondorCollector
fout = io.make_noisy_htcondor_run_fname(0, n, jcoup, dt, tstart, tstop, zeta, eps, n_trajectories)
np.save(os.path.join(dest, fout), out)