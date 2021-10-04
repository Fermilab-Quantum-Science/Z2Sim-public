#include "configgeneration.h"

void generate_configurations(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau, int startconfig)
{
    Lattice lat(ns, nt, coupling, false);
    double action = 0.0;
    //for (int i = 0; i < burninsteps; i++)
    //{
    //    std::cout << "\t mag: " << lat.CalcPlaquetteAverage() << "\n";
    //    lat.Update(((double)i) / ((double)burninsteps) * betatau);
    //}
    for (int i = 0; i < burninsteps; i++)
    {
        lat.Update(betatau);
        std::cout << "\t mag: " << lat.CalcPlaquetteAverage() << "\n";
    }
    if (startconfig != 0)
    {
        lat.LoadConfiguration(startconfig - 1, betatau);
    }
    for (int i = 0; i < configurations; i++) {
        std::cout << i << "\n";
        for (int j = 0; j < montecarlosweeps; j++) {
            lat.Update(betatau);
        }
        lat.SaveConfiguration(i + startconfig, betatau);
    }

}

void generate_configurations_and_calc_plaq(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau, int xshift, bool newstart)
{
	std::string filename;
	std::vector<std::vector<std::string> > plaq_ss_filenames;
	std::vector<std::vector<std::vector<std::string> > > plaq_st_filenames;
	plaq_ss_filenames.resize(ns);
	plaq_st_filenames.resize(ns);
	for (int i = 0; i < ns; i++)
	{
		plaq_st_filenames[i].resize(ns);
	}
	for (int x = 0; x < ns; x++)
	{
		for (int y = 0; y < ns; y++)
		{
			filename = "../PlaquetteLoc/";
			std::string str1 = std::to_string(betatau);
			std::string str2 = std::to_string(coupling);
			str1 = str1.replace(str1.find('.'), 1, "_");
			str2 = str2.replace(str2.find('.'), 1, "_");
			filename += "plaq_ss_x=" + std::to_string(x) + "y=" + std::to_string(y) + "beta=" + str1;
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			filename += ".csv";
			plaq_ss_filenames[x].push_back(filename);
			for (int mu = 0; mu < 2; mu++)
			{
				filename = "../PlaquetteLoc/";
				filename += "plaq_st_x=" + std::to_string(x) + "y=" + std::to_string(y) + "mu" + std::to_string(mu) +  "beta=" + str1;
				filename += "anisotropy=" + str2;
				filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
				filename += ".csv";
				plaq_st_filenames[x][y].push_back(filename);
			}
		}
	}
	// check if we are starting from scratch
	Lattice lat(ns, nt, coupling, false);
	if (newstart)
	{
		std::cout << "starting new lattice....\n";
		std::cout << "Burning in the lattice\n";
		for (int i = 0; i < burninsteps; i++)
		{
			lat.Update(betatau);
		}
		std::cout << "lattice burnin complete\n";
	}
	else
	{
		std::cout << "loading lattice\n";
		lat.LoadConfiguration(0, betatau);
	}
	std::cout << "generating configurations...\n";
	for (int i = 0; i < configurations; i++)
	{
		for (int j = 0; j < montecarlosweeps; j++)
		{
			lat.Update(betatau);
		}
		std::cout << "configuration " << i + 1 << " of " << configurations << " generated\n";
		std::cout << "calculating wilson loops...\n";
		lat.CalculatePlaquetteFromCenter(ns, nt, coupling, betatau, plaq_ss_filenames, plaq_st_filenames);
		std::cout << "wilson loops done!\n" << "calculating correlators...\n";
	}
}

void generate_configurations_and_calc(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau, int xshift, bool newstart, bool calc_corr, bool calc_wil)
{
    // filenames for saving stuff
	std::vector<std::vector<std::string> > wilson_ss_filenames;
	std::vector<std::vector<std::string> > wilson_st_filenames;
	std::vector<std::string> correlator_filenames;
	int s_max = (ns - 1) - 2 * xshift;
	int t_max = nt / 2;
	wilson_ss_filenames.resize(s_max);
	wilson_st_filenames.resize(s_max);
	correlator_filenames.resize(t_max);
	std::string filename;
	for (int length = 0; length < s_max; length++)
	{
		wilson_ss_filenames[length].resize(s_max);
		wilson_st_filenames[length].resize(t_max);
		for (int width = 0; width < s_max; width++)
		{
			filename = "../wilsonloops/";
			std::string str1 = std::to_string(betatau);
			std::string str2 = std::to_string(coupling);
			str1 = str1.replace(str1.find('.'), 1, "_");
			str2 = str2.replace(str2.find('.'), 1, "_");
			filename += "wilsonloop_ss_l=" + std::to_string(length) + "w=" + std::to_string(width) + "beta=" + str1;
			
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			filename += ".csv";
			wilson_ss_filenames[length][width] = filename;
		}
		for (int width = 0; width < t_max; width++)
		{
			filename = "../wilsonloops/";
			std::string str1 = std::to_string(betatau);
			std::string str2 = std::to_string(coupling);
			str1 = str1.replace(str1.find('.'), 1, "_");
			str2 = str2.replace(str2.find('.'), 1, "_");
			filename += "wilsonloop_st_l=" + std::to_string(length) + "w=" + std::to_string(width) + "beta=" + str1;
			
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			filename += ".csv";
			wilson_st_filenames[length][width] = filename;
		}
	}
	for (int width = 0; width < t_max; width++)
	{
		filename = "../correlators/";
		std::string str1 = std::to_string(betatau);
		std::string str2 = std::to_string(coupling);
		str1 = str1.replace(str1.find('.'), 1, "_");
		str2 = str2.replace(str2.find('.'), 1, "_");
		filename += "correlators_w=" + std::to_string(width) + "beta=" + str1;
		filename += "anisotropy=" + str2;
		filename += "ns_" + std::to_string(ns);
		filename += "_nt_" + std::to_string(nt);
		filename += ".csv";
		correlator_filenames[width] = filename;
	}
	for (int length = 0; length < t_max; length++)
	{
		std::cout << correlator_filenames[length] << "\n";
	}
    // check if we are starting from scratch
    Lattice lat(ns, nt, coupling, false);
    if (newstart)
    {
        std::cout << "starting new lattice....\n";
        std::cout << "Burning in the lattice\n";
        for (int i = 0; i < burninsteps; i++)
        {
            lat.Update(betatau);
        }
        std::cout << "lattice burnin complete\n";
    }
    else
    {
        std::cout << "loading lattice\n";
        lat.LoadConfiguration(0, betatau);
    }
    std::cout << "generating configurations...\n";
    for (int i = 0; i < configurations; i++)
    {
        for (int j = 0; j < montecarlosweeps; j++)
        {
            lat.Update(betatau);
        }
        std::cout << "configuration " << i + 1 << " of " << configurations << " generated\n";
		std::cout << "calculating wilson loops...\n";
		lat.CalculateWilsonLoopsOnTheFly(ns, nt, coupling, betatau, xshift, wilson_ss_filenames, wilson_st_filenames, correlator_filenames, calc_corr, calc_wil);
		std::cout << "wilson loops done!\n" << "calculating correlators...\n";
    }
}

void generate_configurations_and_calc_corr(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau)
{

	std::string filename;
	std::vector<std::vector<std::string> > operator_filenames;
	std::vector<std::vector<std::vector<std::string> > > correlator_filenames;
	operator_filenames.resize(28);
	correlator_filenames.resize(28);
	for (int i = 0; i < 28; i++)
	{
		correlator_filenames[i].resize(nt / 2);
		for (int dt = 0; dt < nt / 2; dt++)
		{
			correlator_filenames[i][dt].resize((ns - 1) / 2);
			for (int d = 0; d < (ns - 1) / 2; d++)
			{
				filename = "../correlators/";
				std::string str1 = std::to_string(betatau);
				std::string str2 = std::to_string(coupling);
				str1 = str1.replace(str1.find('.'), 1, "_");
				str2 = str2.replace(str2.find('.'), 1, "_");
				filename += "correlator=" + std::to_string(i);
				filename += "rim=" + std::to_string(d) + "dt=" + std::to_string(dt) + "beta=" + str1;
				filename += "anisotropy=" + str2;
				filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
				filename += ".csv";
				correlator_filenames[i][dt][d] = filename;
			}
		}
		operator_filenames[i].resize((ns - 1) / 2);
		for (int d = 0; d < (ns - 1) / 2; d++)
		{
			filename = "../crosscorrelators/";
			std::string str1 = std::to_string(betatau);
			std::string str2 = std::to_string(coupling);
			str1 = str1.replace(str1.find('.'), 1, "_");
			str2 = str2.replace(str2.find('.'), 1, "_");
			filename += "operator=" + std::to_string(i);
			filename += "rim=" + std::to_string(d) + "beta=" + str1;
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			filename += ".csv";
			operator_filenames[i][d] = filename;
		}
	}// check if we are starting from scratch
	Lattice lat(ns, nt, coupling, false);
	std::cout << "starting new lattice....\n";
	std::cout << "Burning in the lattice\n";
	for (int i = 0; i < burninsteps; i++)
	{
		lat.Update(betatau);
	}
	std::cout << "lattice burnin complete\n";
	std::cout << "generating configurations...\n";
	for (int i = 0; i < configurations; i++)
	{
		for (int j = 0; j < montecarlosweeps; j++)
		{
			lat.Update(betatau);
		}
		std::cout << "configuration " << i + 1 << " of " << configurations << " generated\n";
		std::cout << "calculating wilson loops...\n";
		lat.CalculateCorrelatorsOnTheFly(ns, nt, coupling, betatau, operator_filenames, correlator_filenames);
		std::cout << "wilson loops done!\n" << "calculating correlators...\n";
	}
}

void generate_configurations_and_calc_spacial_wilson_and_corr(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau)
{

	Lattice lat(ns, nt, coupling, false);
	std::cout << "starting new lattice....\n";
	std::cout << "Burning in the lattice\n";
	for (int i = 0; i < burninsteps; i++)
	{
		lat.Update(betatau);
	}
	std::cout << "lattice burnin complete\n";
	std::cout << "generating configurations...\n";
	lat.CalculateWilsonLoopsLocalandCrossCorrelatorsOnTheFly(ns, nt, coupling, betatau, montecarlosweeps, configurations);
}

void generate_configurations_and_calc_spacial_wilson(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau)
{

	Lattice lat(ns, nt, coupling, false);
	std::cout << "starting new lattice....\n";
	std::cout << "Burning in the lattice\n";
	for (int i = 0; i < burninsteps; i++)
	{
		lat.Update(betatau);
	}
	std::cout << "lattice burnin complete\n";
	std::cout << "generating configurations...\n";
	lat.CalculateWilsonLoopsLocalOnTheFly(ns, nt, coupling, betatau, montecarlosweeps, configurations);
}

void generate_configurations_and_calc_cross_corr(int ns, int nt, int burninsteps, int montecarlosweeps, int configurations, double coupling, double betatau)
{
	std::string filename;
	std::vector<std::vector<std::string> > operator_filenames;
	std::vector<std::vector<std::vector<std::vector<std::string> > > > correlator_filenames;
	operator_filenames.resize(28);
	correlator_filenames.resize(28);
	for (int i = 0; i < 28; i++)
	{
		correlator_filenames[i].resize(28);
		for (int j = 0; j < 28; j++)
		{
			correlator_filenames[i][j].resize(nt / 2);
			for (int dt = 0; dt < nt / 2; dt++)
			{
				correlator_filenames[i][j][dt].resize((ns - 1) / 2);
				for (int d = 0; d < (ns - 1) / 2; d++)
				{
					filename = "../crosscorrelators/";
					std::string str1 = std::to_string(betatau);
					std::string str2 = std::to_string(coupling);
					str1 = str1.replace(str1.find('.'), 1, "_");
					str2 = str2.replace(str2.find('.'), 1, "_");
					filename += "correlator=" + std::to_string(i) + "_" + std::to_string(j);
					filename += "rim=" + std::to_string(d) + "dt=" + std::to_string(dt) + "beta=" + str1;
					filename += "anisotropy=" + str2;
					filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
					filename += ".csv";
					correlator_filenames[i][j][dt][d] = filename;
				}
			}
		}
		operator_filenames[i].resize((ns - 1) / 2);
		for (int d = 0; d < (ns - 1) / 2; d++)
		{
			filename = "../crosscorrelators/";
			std::string str1 = std::to_string(betatau);
			std::string str2 = std::to_string(coupling);
			str1 = str1.replace(str1.find('.'), 1, "_");
			str2 = str2.replace(str2.find('.'), 1, "_");
			filename += "operator=" + std::to_string(i);
			filename += "rim=" + std::to_string(d) + "beta=" + str1;
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			filename += ".csv";
			operator_filenames[i][d] = filename;
		}
	}// check if we are starting from scratch
	Lattice lat(ns, nt, coupling, false);
	std::cout << "starting new lattice....\n";
	std::cout << "Burning in the lattice\n";
	for (int i = 0; i < burninsteps; i++)
	{
		lat.Update(betatau);
	}
	std::cout << "lattice burnin complete\n";
	std::cout << "generating configurations...\n";
	for (int i = 0; i < configurations; i++)
	{
		for (int j = 0; j < montecarlosweeps; j++)
		{
			lat.Update(betatau);
		}
		std::cout << "configuration " << i + 1 << " of " << configurations << " generated\n";
		std::cout << "calculating wilson loops...\n";
		lat.CalculateCrossCorrelatorsOnTheFly(ns, nt, coupling, betatau, operator_filenames, correlator_filenames);
		std::cout << "wilson loops done!\n" << "calculating correlators...\n";
	}
}
