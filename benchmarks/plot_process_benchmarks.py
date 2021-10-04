import numpy as np
import argparse

import matplotlib.pyplot as plt

"""
The following data are taken from some diagnostic runs on t100 gpus with
noise trajectories.
"""

# Fix noise and physical parameters
#
tags = [
    (0, "c2_standard_30"),
    (0, "n1_standard_32"),
    (1, "gpu_tesla_a100"),
    (1, "gpu_tesla_p100"),
    (1, "gpu_tesla_t4"),
    (1, "gpu_tesla_v100")
]
# Nested dictionary of grid sizes vs. trajectory counts
fig, axes = plt.subplots(2, 3, sharey=True, sharex=True)
x0, x1 = 1, 100

grid_sweep = [3, 4, 5]

colors = ['r', 'g', 'b']

for j, (gpu, fstr) in enumerate(tags):
    # 20 trotter steps for x trajectories
    results = np.load(f"./benchmark_results/benchmark_gpu{gpu}_{fstr}_tstart20.npy")
    ax = axes.flatten()[j]
    for i, n in enumerate(grid_sweep):
        y0 = results[0,i]
        y1 = results[1,i]
        # *5 because this was 20 trottersteps in the simulation
        slope = 5 * (y1 - y0)  /  (x1 - x0)

        label = f'{n}x{n}; m={slope:3.5f} s/(100 t*ts); b={y0:3.2f} s'
        ls = '--' if gpu else '-'
        ax.plot([x0, x1], [y0, y1], ls=ls, label=label, c=colors[i])

    ax.set_title(fstr, size=20)

    ax.set_xlabel("number of trajectories")
    ax.set_ylabel("wall time (s)")
    ax.set_xlim(x0, x1)
    ax.set_ylim(0, ax.get_ylim()[1])
    # ax.set_title(f"GPU={gpu}", size=22)
    ax.legend(loc="upper left")
plt.show()