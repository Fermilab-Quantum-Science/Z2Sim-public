
#ifndef WILSONHEADERDEF
#define WILSONHEADERDEF

#include "Lattice.h"
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include <omp.h>

void CalculateWilsonLoopsFromConfig(int ns, int nt, int configstart, int number_of_configs, double coupling, double betatau, int set);

void CalculateWilsonLoopsOnTheFly(Lattice lat, int ns, int nt, double coupling, double betatau, int xshift, std::vector<std::vector< std::string> > wilson_ss_filenames, std::vector<std::vector< std::string> > wilson_st_filenames, std::vector<std::string> correlator_filenames);
#endif