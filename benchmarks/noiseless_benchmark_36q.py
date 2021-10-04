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
parser.add_argument('-resource', metavar='resource', type=str, nargs=1,
                    help='string tag describing this device')
args = parser.parse_args()

resource = args.resource[0]
TROTTER_STEPS = 10
jcoup = 0.8
dt = 0.25

# Initialize simulators and noise model
N_FUSE = 4
# # # # # # # # # # # # # # # # # # # # # # # # #
# CHECK THIS BEFOR RUNNING LOCALLY
N_THREADS = 160
n = 6
# # # # # # # # # # # # # # # # # # # # # # # # #

qsim_options = qsimcirq.QSimOptions(
    max_fused_gate_size = N_FUSE,
    cpu_threads = N_THREADS,
    denormals_are_zeros = True,
)

# "SIM A" ancillafree observable computation using 4x 36-qubit circuit
# The outermost loop in this routine is over the four distinct circuits
# So the first 25% of debug statement timing results will correspond
# to a single choice of initial generator
out, times_out = production.compute_ancillafree_obs_with_intermediate_state_vector(
    n=n,
    trotter_steps=TROTTER_STEPS,
    jcoup=jcoup,
    dt=dt,
    qsim_options=qsim_options,
    all_observables=True,
    obc=True,
    decompbasis="MS",
    verbose=True,
)

# OLD simulation using extra ancilla qubit.
# This has been updated to only compute 21/36 observables.
# out, times_out = production.compute_obs_with_intermediate_state_vector(
#     n=n,
#     trotter_steps=TROTTER_STEPS,
#     jcoup=jcoup,
#     dt=dt,
#     qsim_options=qsim_options,
#     all_observables=True,
#     obc=True,
#     decompbasis="MS",
#     verbose=True,
# )

np.save(f"benchmark_{resource}_trotter{TROTTER_STEPS}.npy", times_out)