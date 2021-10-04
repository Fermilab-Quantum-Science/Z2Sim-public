#pragma once
#include <iostream>
#include <stdlib.h>
#include "Lattice.h"

void MeasureCorrelatorFromConfig(int ns, int nt, int configstart, int number_of_configs, double coupling, double betatau, int set);

void MeasureCorrelatorOnTheFly(Lattice lat, int ns, int nt, double coupling, double betatau, int xshift, std::vector<std::string> correlator_filenames);
