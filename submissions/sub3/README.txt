Production level run for n=3

I've sliced the rectangular noise parameter sweep into 3 subrectangles. Each
of these will be submitted as a separate job, partly in case anything goes
wrong during a single long job, and partly to parallelize across more gpu
resources.