
# Z2Sim

This repository is for the Google Cloud Burst and Fermilab collaboration on large lattice simulations of 2 + 1-D **Z** <sub>2</sub>
gauge theory using both dual and gauge representations of the theory.

## Installation

The package `z2_sim` should be installed as a local package in order to run the scripts and tests:
```
cd Z2Sim-public
python3 -m pip install -e z2_sim
```
Note that this will install requirements automatically. However, any GPU usage will require `qsimcirq` to be installed from source and appropriate graphics card drivers.

## Noise model usage

Noisy operations are implemented via cirq's `cirq.devices.NoiseModel` interface. We have used several noise models based on uniform error across all qubits. To initialize a noise model with uniform errors, just define a parameter describing typical error and then initialize the desired noise model. For example, we used a crosstalk model of the form:
```
from z2_sim.src.QuantumCircuits.Cirq_Code.noise.zz_crosstalk_v1 import ZZCrossTalkV1
# Initialize a noise model that applies a ZZ rotation of the form
# exp( (-i2π ζT)* |11><11|  ) after every sqrtISWAP gate
zeta = 150000 # ζ
target_gate = "SIS" # T
GATE_DURATION = 1e-8
zeta_model = ZZCrossTalkV1(zeta, target_gate, gate_duration=GATE_DURATION, sampled_error=False)
```
We composed this with a two-local depolarization model of the form
```
from z2_sim.src.QuantumCircuits.Cirq_Code.noise.two_local_depol import TwoLocalDepol
# Initialize a noise model with one- and two-local symmetric depolarizing noise
# with the depolarization rate on two-qubit gates being 10x stronger than on
# one-qubit gates
eps_noise_model = TwoLocalDepol(err_1q=eps / 10, err_2q=eps)
```

Note that these noise models are designed to be composable using cirq's tagging feature. You should always inspect composed noise models to make sure that noisy operations are not being applied as errors in another noise operation, e.g. we do not want the second noise model to apply two-local depolarization to any ZZ gate introduced by the first noise model.

Once noise models are initialized, they can be applied to an existing circuit like this:
```
circuit = cirq.Circuit(...)
for noise_model in [zeta_model, eps_noise_model]:
    noisy_circuit = circuit.with_noise(noise_model)
```

Note that in general the order of composition matters if the noisy operations are noncommuting, but in the case of depolarizing noise it does not matter whether it is applied before or after a gate in an isolated moment.

Once a noisy circuit is constructed, it can either be simulated by cirq's built-in desity matrix simulator `cirq.DensityMatrixSimulator` or by performing trajectory simulations on qsim:
```
# qsim trajectory simulation (approximate)
n_trajectories = 10_000 # number of noisy trajectories
observables = [...]

qsim_options = qsimcirq.QSimOptions(
    ev_noisy_repetitions=n_trajectories,
)
qsim_simulator = qsimcirq.QSimSimulator(qsim_options=qsim_options)

expectation_values = qsim_simulator.simulate_expectation_values(
    noisy_circuit,
    observables=observables,
)

# cirq density matrix simulator (floating point precision)
# This will be exponentially more expensive than trajectory simulation in
# general but may be practical for small systems where the trajectory count
# dominates simulation time.
cirq_simulator = cirq.DensityMatrixSimulator()

expectation_values = cirq_simulator..simulate_expectation_values(
    noisy_circuit,
    observables=[XX_obs]
)
```
