Production level run for n=4 (and maybe n=3 repeats if desired)

The main difference here is that the condor jobs will be distributed
with _different parameter values_. So instead of duplicating a single 100-shot
sweep 10 times, I distribute many 1000-shot sweeps involving different
noise strengths.

This is motivated by the large overhead for GPU simulations. This different
submission style is accomplished by a hack wherein a shell script accepts the
condor proc id as an input, and then dispatches a python exe with a different
parameter set depending on that ID.

The condor outputs can still be collected in the same way, since I will hardcode
the output filename to have proc=0. There will be no name clashes since the
individual jobs will output a unique file signature based on the hardcoded
phyiscal and noise parameters.