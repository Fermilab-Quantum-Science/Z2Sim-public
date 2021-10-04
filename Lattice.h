#ifndef Z2HEADERDEF
#define Z2HEADERDEF

#include <math.h>
#include <time.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

class Lattice
{
public:
	// actual lattice
	short int* link;
	// array for shuffling index elements
	unsigned int* indexarray;
	// wilson loop averages for config
	double* wilson_ss;
	double* wilson_st;
	// wilson loop norms
	double* norms_ss;
	double* norms_st;
	/*==========================================*/

	// class constructor
	Lattice(int spacial_dim, int temporal_dim, double coup, bool randomstart);
	// class destructor
	~Lattice();

	// move the vector x 1 step forward in direction d
	void MoveUp(int x[], int d);

	// move the vector x, 1 step backward in direction d
	void MoveDown(int x[], int d);

	// get the index of the array for the given vector
	int GetIndex(int x[], int d);

	// shuffle the elements of the indexarray
	void Shuffle();

	// monte carlo update across all sites at temporal temperature betat
	void Update(double betat);

	// calculate the average plaquette value
	double CalcPlaquetteAverage();

	// function to save the configuration
	void SaveConfiguration(int configuration_number, double betat);

	// function to load the configuration
	void LoadConfiguration(int configuration_number, double betat);

	// reset wilson loop
	void ResetWilsonLoops();

	// generate norms
	void CalculateNorms();


	// calculate the cross correlator terms
	std::vector<std::vector<double> > CalculateOperators(int t);

	void CalculateCrossCorrelatorsOnTheFly(int ns, int nt, double coupling, double betatau, std::vector<std::vector<std::string> > operator_filenames,
		std::vector<std::vector<std::vector<std::vector<std::string> > > > correlatorfilenames);

	void CalculateWilsonLoopsLocalandCrossCorrelatorsOnTheFly(int ns, int nt, double coupling, double betatau, int montecarlosweeps, int configurations);

	void CalculateWilsonLoopsLocalOnTheFly(int ns, int nt, double coupling, double betatau, int montecarlosweeps, int configurations);

	void CalculateCorrelatorsOnTheFly(int ns, int nt, double coupling, double betatau, std::vector<std::vector<std::string> > operator_filenames,
		std::vector<std::vector<std::vector<std::string> > > correlatorfilenames);

	// calculate a height x width wilson loop
	int CalculateWilsonLoop(int x[], int d, int dperp, int height, int width);

	// calculate all available spacial wilson loops
	void CalculateAllSpacialLoops();

	// calculate all the available temporal loops
	void CalculateAllTemporalLoops();

	void CalculateWilsonLoopsOnTheFly(int ns, int nt, double coupling, double betatau, int xshift, std::vector<std::vector<std::string>> wilson_ss_filenames, std::vector<std::vector<std::string>> wilson_st_filenames, std::vector<std::string> correlator_filenames,
		                              bool calc_corr, bool calc_wil);

	void CalculatePlaquetteFromCenter(int ns, int nt, double coupling, double betatau, std::vector<std::vector<std::string> > plaqss_loc_filenames, std::vector<std::vector<std::vector<std::string> > >plaqst_loc_filenames);

	// calculate Wilson Loops on the fly
	void CalculateWilsonLoops(int confignum, double betatau);

private:
	// spacial (ns) and (nt) temporal length of the lattice
	int ns, nt;
	int dim;
	// magnetic field coupling
	double coupling;
};

#endif


