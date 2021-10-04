
# Z2Sim

This repository is for the Google Cloud Burst and Fermilab collaboration on large lattice simulations of 2 + 1-D **Z** <sub>2</sub>
gauge theory using both dual and gauge representations of the theory.

## Installation notes

Please note that this repository has been structured with specific applications in mind. The `requirements.txt` is meant for use in constructing the Dockerfile based off the latest qsim build at gcr.io. Therefore a local installation of this repository may not work properly.

## Overview

The Hamiltonian describing using the gauge representation is:

![](docs/imgs/hgauge.png)

The Hamiltonian for the dual representation is,

![](docs/imgs/hdual.png)

The **Z**<sub>2</sub> gauge theory has two quantum phases: a confined (J >> &Gamma;) and deconfined (&Gamma; >> J) phase.


None of the layout here is set in stone and just an attempt to make some organization.


## Repository Layout

The source code is laid out into several folders:
1. epoq-Monte-Carlo
	- This has the code to do the euclidean monte carlo simulations to generate samples from the denisty  matrix to
2. NumericsPython
	- This has the numpy code to carry out sparse matrix evolution which can be used as a consistency check for the cirq code
3. QuantumCircuit
	- This folder contains the code for both cirq and QISKit code to implement the quantum circuits for the various lattices


## Current To-Do list

Several things need to currently be implemented with respect to coding for the cirq classes, numeric simulation code to carry out Trotterization noiselessly so that we have comparisons for small lattice simulations, and possibly other avenues if we are interested in certain types of state preparation.


In regard to physics we need to look at the quantum phases and find the set of parameters that give us the continuum limit.

The to do list is not comprehensive and just what needs to be done with (Erik's) code


-------------

### Coding

-----------


#### Cirq code

-  [x] ~~Optimize the CNOT gates into either `FSIM` or `sqrt-iswap` gates~~
-  [ ] Develop circuit decompositions for the Moeller  in terms of arbitrary
-  [ ] Develop algorithm to generate initial states for the dual representations
-  [ ] Develop algorithm to generate initial states for the gauge representation
-  [ ] Develop algorithm to force gauge invariant initial state
-  [ ] Develop method for state preparation
-  [ ]


#### Noiseless Code

- [ ] Add code to simulate arbitrary initial states

-----------------------

### Physics

-----------------------

#### Simulation

- [ ] What are the dynamics / degrees of freedom in each phase?
- [x] ~~What dynamical quantities are we interested in?~~ We are going to look at Wilson loop correlators <0|W_ij(x,t)W_kl(y,0)|0>
- [ ] How do we extract the lattice scales a_s and a_t from the Wilson loop correlators?
- [ ] What kind of state preparation do we want to do? Naive state preparation seems to be the one that is most feasible

#### Correlators

We need to set the scale by using the correlators

