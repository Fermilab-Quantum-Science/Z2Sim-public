#ifndef CONFIGGENHEADERDEF
#define CONFIGGENHEADERDEF

#include "Lattice.h"
#include "WilsonLoopCalc.h"
#include "Correlator.h"

void generate_configurations(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau, int startconfig);


void generate_configurations_and_calc_plaq(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau, int xshift, bool newstart);

void generate_configurations_and_calc(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau, int xshift, bool newstart, bool calc_corr, bool calc_wil);

void generate_configurations_and_calc_cross_corr(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau);

void generate_configurations_and_calc_corr(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau);

void generate_configurations_and_calc_spacial_wilson(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau);

void generate_configurations_and_calc_spacial_wilson_and_corr(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau);
#endif