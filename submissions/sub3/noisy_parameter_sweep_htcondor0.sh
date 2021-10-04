#!/bin/bash
# You must source the installed drivers and libraries as they appear in the
# terraform startup script.
export PATH="/usr/local/cuda-11.4/bin${PATH:+:${PATH}}"
export PYTHONPATH=$PYTHONPATH:"/home/peterse583/qsim"

# Install the wheel that was passed in by htcondor.
# Assumes that the dependencies are already on the machine...
sudo python3 -m pip install z2_sim-0.1-py3-none-any.whl

# If this directory already exists life goes on just fine.
DEST="./temp"
mkdir $DEST

# python3 noisy_parameter_sweep_htcondor.py -proc $1 -n 3 -j 0.714285 0.625 .555556 -dt .25 -tstart 1 -tstop 51 -zeta 0 -eps 0.001 0.0001 0.00001 -r 100 -dest $DEST
python3 noisy_parameter_sweep_htcondor.py -proc $1 -n 3 -j 0.714285 0.625 .555556 \
  -dt .25 -tstart 1 -tstop 51 \
  -zeta 300000 450000 \
  -eps 0.001 0.002 0.003 0.004 0.005 -r 100 -dest $DEST