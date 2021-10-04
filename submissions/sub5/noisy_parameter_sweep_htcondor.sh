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

# Modify the parameters within the below file.
python3 noisy_parameter_sweep_htcondor.py -proc $1 -dest $DEST
