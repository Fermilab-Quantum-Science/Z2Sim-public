#include "Correlator.h"

void MeasureCorrelatorFromConfig(int ns, int nt, int configstart, int number_of_configs, double coupling, double betatau, int set)
{
	int xshift = (ns - 1) / 4;
	int x[3];
	std::vector<std::string> correlator_filenames;
	correlator_filenames.resize(nt / 2);
	std::string filename;
	for (int width = 0; width < nt / 2; width++)
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
	Lattice lat(ns, nt, coupling, false);
	for (int configuration_num = configstart; configuration_num < number_of_configs; configuration_num++)
	{
		std::cout << "starting configuration " << configuration_num << "...\n";
		std::vector<std::vector<int> > correlators;
		std::vector<double> means;
		std::vector<double> sdevs;
		means.resize(nt / 2);
		sdevs.resize(nt / 2);
		correlators.resize(nt / 2);
		lat.LoadConfiguration(configuration_num, betatau);
		//std::cout << "configuration loaded\n";
		// iterate through the x sites away from 1/4 of the edges x side
		for (x[0] = xshift; x[0] < xshift + (ns - 1) / 2; x[0]++)
		{
			// iterate through the y sites away from 1/4 of the edges y side
			for (x[1] = xshift; x[1] < xshift + (ns - 1) / 2; x[1]++)
			{
				// iterate through the time steps
				for (x[2] = 0; x[2] < nt; x[2]++)
				{
					for (int w = 0; w < nt / 2; w++)
					{
						int y[3];
						int z[3];
						y[0] = 0 + x[0];
						y[1] = 0 + x[1];
						y[2] = 0 + x[2];
						z[0] = 0 + x[0];
						z[1] = 0 + x[1];
						z[2] = (w + x[2]) % nt;
						int t1 = lat.CalculateWilsonLoop(y, 0, 1, 1, 1);
						int t2 = lat.CalculateWilsonLoop(z, 0, 1, 1, 1);
						correlators[w].push_back(t1 * t2);
					}
				}
			}
		}
		std::cout << "correlators calculated\n";
		double val = 0.0;
		for (int w = 0; w < nt / 2; w++)
		{
			double num_elements = ((double)correlators[w].size());
			for (int i = 0; i < num_elements; i++)
			{
				val = ((double)correlators[w][i]);
				if (val > 1.0)
				{
					val = 1.0;
				}
				if (val < -1.0)
				{
					val = -1.0;
				}
				means[w] += val;
			}
			means[w] /= num_elements;
			double variance = 0.0;
			for (int i = 0; i < num_elements; i++)
			{
				val = ((double)correlators[w][i]);
				if (val > 1.0)
				{
					val = 1.0;
				}
				if (val < -1.0)
				{
					val = -1.0;
				}
				variance += (val - means[w]) * (val - means[w]);
			}
			sdevs[w] = sqrt(variance / num_elements);
		}
		std::cout << "Saving data....";
		std::ofstream outfile;
		int counter = 0;
		for (int w = 0; w < nt / 2; w++)
		{
			outfile.open(correlator_filenames[w], std::ios::app);
			outfile << means[w] << "," << sdevs[w] << "\n";
			outfile.close();
		}
		std::cout << "configuration done!\n";
	}
}

void MeasureCorrelatorOnTheFly(Lattice lat, int ns, int nt, double coupling, double betatau, int xshift, std::vector<std::string> correlator_filenames)
{
	int x[3];
	std::vector<std::vector<int> > correlators;
	std::vector<double> means;
	std::vector<double> sdevs;
	means.resize(nt / 2);
	sdevs.resize(nt / 2);
	correlators.resize(nt / 2);
	// iterate through the x sites away from 1/4 of the edges x side
	for (x[0] = xshift; x[0] < xshift + (ns - 1) / 2; x[0]++)
	{
		// iterate through the y sites away from 1/4 of the edges y side
		for (x[1] = xshift; x[1] < xshift + (ns - 1) / 2; x[1]++)
		{
			// iterate through the time steps
			for (x[2] = 0; x[2] < nt; x[2]++)
			{
				for (int w = 0; w < nt / 2; w++)
				{
					std::cout << x[0] << "\n";
					int y[3];
					int z[3];
					y[0] = 0 + x[0];
					y[1] = 0 + x[1];
					y[2] = 0 + x[2];
					z[0] = 0 + x[0];
					z[1] = 0 + x[1];
					z[2] = (w + x[2]) % nt;
					int t1 = lat.CalculateWilsonLoop(x, 0, 1, 1, 1);
					int t2 = lat.CalculateWilsonLoop(z, 0, 1, 1, 1);
					//correlators[w].push_back(t1 * t2);
					std::cout << t1  << t2 << "\n";
				}
			}
		}
	}
	std::cout << "correlators calculated\n";
	double val = 0.0;
	for (int w = 0; w < nt / 2; w++)
	{
		double num_elements = ((double)correlators[w].size());
		for (int i = 0; i < num_elements; i++)
		{
			val = ((double)correlators[w][i]);
			if (val > 1.0)
			{
				val = 1.0;
			}
			if (val < -1.0)
			{
				val = -1.0;
			}
			means[w] += val;
		}
		means[w] /= num_elements;
		double variance = 0.0;
		for (int i = 0; i < num_elements; i++)
		{
			val = ((double)correlators[w][i]);
			if (val > 1.0)
			{
				val = 1.0;
			}
			if (val < -1.0)
			{
				val = -1.0;
			}
			variance += (val - means[w]) * (val - means[w]);
		}
		sdevs[w] = sqrt(variance / num_elements);
	}
	std::cout << "Saving data....";
	std::ofstream outfile;
	int counter = 0;
	for (int w = 0; w < nt / 2; w++)
	{
		std::cout << correlator_filenames[w] << "\n";
		outfile.open(correlator_filenames[w], std::ios::app);
		outfile << means[w] << "," << sdevs[w] << "\n";
		outfile.close();
	}
}
