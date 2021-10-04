#include "WilsonLoopCalc.h"



void CalculateWilsonLoopsFromConfig(int ns, int nt, int configstart, int number_of_configs, double coupling, double betatau, int set)
{
	// construct the lattice object
	int xshift = (ns - 1) / 4;
	int x[3];
	std::vector<std::vector<std::string> > wilson_ss_filenames;
	std::vector<std::vector<std::string> > wilson_st_filenames;
	wilson_ss_filenames.resize(ns / 2);
	wilson_st_filenames.resize(nt / 2);
	std::string filename;
	for (int length = 0; length < ns / 2; length++)
	{
		for (int width = 0; width < ns / 2; width++)
		{
			filename = "../wilsonloops/";
			std::string str1 = std::to_string(betatau);
			std::string str2 = std::to_string(coupling);
			str1 = str1.replace(str1.find('.'), 1, "_");
			str2 = str2.replace(str2.find('.'), 1, "_");
			if (width > length)
			{
				filename += "wilsonloop_ss_l=" + std::to_string(length) + "w=" + std::to_string(width) + "beta=" + str1;
			}
			else
			{
				filename += "wilsonloop_ss_l=" + std::to_string(width) + "w=" + std::to_string(length) + "beta=" + str1;
			}
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			if (set != 0)
			{
				filename += "set=" + std::to_string(set);
			}
			filename += ".csv";
			wilson_ss_filenames[length].push_back(filename);
		}
	}
	for (int length = 0; length < nt / 2; length++)
	{
		for (int width = 0; width < ns / 2; width++)
		{
			filename = "../wilsonloops/";
			std::string str1 = std::to_string(betatau);
			std::string str2 = std::to_string(coupling);
			str1 = str1.replace(str1.find('.'), 1, "_");
			str2 = str2.replace(str2.find('.'), 1, "_");
			filename += "wilsonloop_st_l=" + std::to_string(width) + "w=" + std::to_string(length) + "beta=" + str1;
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			if (set != 0)
			{
				filename += "set=" + std::to_string(set);
			}
			filename += ".csv";
			wilson_st_filenames[length].push_back(filename);
		}
	}
	Lattice lat(ns, nt, coupling, false);
	for (int configuration_num = configstart; configuration_num < number_of_configs; configuration_num++)
	{
		std::cout << "starting configuration " << configuration_num << "...";
		// construct an object to hold the spacial-spacial wilson loops values
		std::vector<std::vector<std::vector<int> > > loops_ss;
		std::vector<std::vector<double> > means_ss;
		std::vector<std::vector<double> > sigmas_ss;
		// construct an object to hold the spacial-temporal wilson loops values;
		std::vector<std::vector<std::vector<int> > > loops_st;
		std::vector<std::vector<double> > means_st;
		std::vector<std::vector<double> > sigmas_st;

		// adjust the vectors to be the right size
		int x[3];
		int l, w;
		loops_ss.resize(ns / 2);
		loops_st.resize(nt / 2);
		means_ss.resize(ns / 2);
		sigmas_ss.resize(ns / 2);
		means_st.resize(nt / 2);
		sigmas_st.resize(nt / 2);
		for (int i = 0; i < ns / 2; i++)
		{
			loops_ss[i].resize(ns / 2);
			means_ss[i].resize(ns / 2);
			sigmas_ss[i].resize(ns / 2);
		}
		for (int i = 0; i < nt / 2; i++)
		{
			loops_st[i].resize(ns / 2);
			means_st[i].resize(ns / 2);
			sigmas_st[i].resize(ns / 2);
		}
		lat.LoadConfiguration(configuration_num, betatau);
		std::cout << "configuration loaded\n";
		int counter1 = 0;
		// iterate through the x sites away from 1/4 of the edges x side
		for (x[0] = xshift; x[0] < xshift + (ns - 1) / 2; x[0]++)
		{
			// iterate through the y sites away from 1/4 of the edges y side
			for (x[1] = xshift; x[1] < xshift + (ns - 1) / 2; x[1]++)
			{
				// iterate through the time steps
				for (x[2] = 0; x[2] < nt; x[2]++)
				{
					// iterate through the possible widthS
					for (int w = 0; w < (ns - 1) / 2 - (x[0] - xshift); w++)
					{
						// iterate through the possible lengths
						for (int l = 0; l < (ns - 1) / 2 - (x[1] - xshift); l++)
						{
							int y[3];
							y[0] = 0 + x[0];
							y[1] = 0 + x[1];
							y[2] = 0 + x[2];
							if (l > w)
							{
								loops_ss[w][l].push_back(lat.CalculateWilsonLoop(y, 0, 1, l + 1, w + 1));
							}
							else
							{
								loops_ss[l][w].push_back(lat.CalculateWilsonLoop(y, 0, 1, l + 1, w + 1));
							}
						}
					}
					for (int l = 0; l < (ns - 1) / 2 - (x[0] - xshift); l++)
					{
						for (int w = 0; w < nt / 2; w++)
						{
							int y[3];
							y[0] = 0 + x[0];
							y[1] = 0 + x[1];
							y[2] = 0 + x[2];
							loops_st[w][l].push_back(lat.CalculateWilsonLoop(y, 0, 2, l + 1, w + 1));
						}
					}
					for (int l = 0; l < (ns - 1) / 2 - (x[1] - xshift); l++)
					{
						for (int w = 0; w < nt / 2; w++)
						{
							int y[3];
							y[0] = 0 + x[0];
							y[1] = 0 + x[1];
							y[2] = 0 + x[2];
							loops_st[w][l].push_back(lat.CalculateWilsonLoop(y, 1, 2, l + 1, w + 1));
						}
					}
				}
			}
		}
		std::cout << "loops calculated\n";
		double val = 0.0;
		for (int w = 0; w < (ns - 1) / 2; w++)
		{
			for (int l = 0; l <= w; l++)
			{
				double num_elements = ((double)loops_ss[l][w].size());
				// calculate the mean
				for (int i = 0; i < num_elements; i++)
				{
					val = ((double)loops_ss[l][w][i]);
					if (val > 1.0)
					{
						val = 1.0;
					}
					else if (val < -1.0)
					{
						val = -1.0;
					}
					means_ss[l][w] += val;
				}
				means_ss[l][w] /= num_elements;
				// calculate the uncertainty
				double variance = 0.0;
				for (int i = 0; i < num_elements; i++)
				{
					val = ((double)loops_ss[l][w][i]);
					if (val > 1.0)
					{
						val = 1.0;
					}
					else if (val < -1.0)
					{
						val = -1.0;
					}
					variance += (val - means_ss[l][w]) * (val - means_ss[l][w]) / num_elements;
				}
				sigmas_ss[l][w] = sqrt(variance / num_elements);
			}
		}
		for (int w = 0; w < (ns - 1) / 2; w++)
		{
			for (int l = 0; l < nt / 2; l++)
			{
				double num_elements = ((double)loops_st[l][w].size());
				// calculate the mean
				for (int i = 0; i < num_elements; i++)
				{
					val = ((double)loops_st[l][w][i]);
					if (val > 1.0)
					{
						val = 1.0;
					}
					else if (val < -1.0)
					{
						val = -1.0;
					}
					means_st[l][w] += val;
				}
				means_st[l][w] /= num_elements;
				// calculate the uncertainty
				double variance = 0.0;
				for (int i = 0; i < num_elements; i++)
				{
					val = ((double)loops_st[l][w][i]);
					if (val > 1.0)
					{
						val = 1.0;
					}
					else if (val < -1.0)
					{
						val = -1.0;
					}
					variance += (val - means_st[l][w]) * (val - means_st[l][w]) / num_elements;
				}
				sigmas_st[l][w] = sqrt(variance / num_elements);
			}
		}
		/*for (int w = 0; w < (ns - 1) / 2; w++)
		{
			for (int l = 0; l <= w; l++)
			{
				means_all_ss[l + (ns - 1) / 2 * w].push_back(means_ss[l + (ns - 1) / 2 * w]);
				sdevs_all_ss[l + (ns - 1) / 2 * w].push_back(sigmas_ss[l + (ns - 1) / 2 * w]);
			}
			for (int l = 0; l < nt / 2; l++)
			{
				means_all_st[l + nt / 2 * w].push_back(means_st[l + nt / 2 * w]);
				sdevs_all_st[l + nt / 2 * w].push_back(sigmas_st[l + nt / 2 * w]);
			}
		}*/
		std::cout << "saving data\n";
		std::ofstream outfile;
		int counter = 0;
		for (int w = 0; w < (ns - 1) / 2; w++)
		{
			for (int l = 0; l <= w; l++)
			{
				outfile.open(wilson_ss_filenames[l][w], std::ios::app);
				outfile << means_ss[l][w] << "," << sigmas_ss[l][w] << "\n";
				outfile.close();
			}
			for (int l = 0; l < nt / 2; l++)
			{
				//counter += 1;
				outfile.open(wilson_st_filenames[l][w], std::ios::app);
				outfile << means_st[l][w] << "," << sigmas_st[l][w] << "\n";
				outfile.close();
			}
		}
		std::cout << "done!\n";
	}
}

void CalculateWilsonLoopsOnTheFly(Lattice lat, int ns, int nt, double coupling, double betatau, int xshift, std::vector<std::vector< std::string> > wilson_ss_filenames, std::vector<std::vector< std::string> > wilson_st_filenames, std::vector<std::string> correlator_filenames)
{
	if (1 == 1)
	{
		int x[3];
		std::vector<std::vector<std::vector<int> > > loops_ss;
		std::vector<std::vector<double> > means_ss;
		std::vector<std::vector<double> > sigmas_ss;
		// construct an object to hold the spacial-temporal wilson loops values;
		std::vector<std::vector<std::vector<int> > > loops_st;
		std::vector<std::vector<double> > means_st;
		std::vector<std::vector<double> > sigmas_st;

		// adjust the vectors to be the right size
		int l, w;
		loops_ss.resize(ns / 2);
		loops_st.resize(nt / 2);
		means_ss.resize(ns / 2);
		sigmas_ss.resize(ns / 2);
		means_st.resize(nt / 2);
		sigmas_st.resize(nt / 2);
		for (int i = 0; i < ns / 2; i++)
		{
			loops_ss[i].resize(ns / 2);
			means_ss[i].resize(ns / 2);
			sigmas_ss[i].resize(ns / 2);
		}
		for (int i = 0; i < nt / 2; i++)
		{
			loops_st[i].resize(ns / 2);
			means_st[i].resize(ns / 2);
			sigmas_st[i].resize(ns / 2);
		}
		std::cout << "configuration loaded\n";
		int counter1 = 0;
		// iterate through the x sites away from 1/4 of the edges x side
		for (x[0] = xshift; x[0] < xshift + (ns - 1) / 2; x[0]++)
		{
			// iterate through the y sites away from 1/4 of the edges y side
			for (x[1] = xshift; x[1] < xshift + (ns - 1) / 2; x[1]++)
			{
				// iterate through the time steps
				for (x[2] = 0; x[2] < nt; x[2]++)
				{
					// iterate through the possible widthS
					for (int w = 0; w < (ns - 1) / 2 - (x[0] - xshift); w++)
					{
						// iterate through the possible lengths
						for (int l = 0; l < (ns - 1) / 2 - (x[1] - xshift); l++)
						{
							int y[3];
							y[0] = 0 + x[0];
							y[1] = 0 + x[1];
							y[2] = 0 + x[2];
							if (l > w)
							{
								loops_ss[w][l].push_back(lat.CalculateWilsonLoop(y, 0, 1, l + 1, w + 1));
							}
							else
							{
								loops_ss[l][w].push_back(lat.CalculateWilsonLoop(y, 0, 1, l + 1, w + 1));
							}
						}
					}
					for (int l = 0; l < (ns - 1) / 2 - (x[0] - xshift); l++)
					{
						for (int w = 0; w < nt / 2; w++)
						{
							int y[3];
							y[0] = 0 + x[0];
							y[1] = 0 + x[1];
							y[2] = 0 + x[2];
							loops_st[w][l].push_back(lat.CalculateWilsonLoop(y, 0, 2, l + 1, w + 1));
						}
					}
					for (int l = 0; l < (ns - 1) / 2 - (x[1] - xshift); l++)
					{
						for (int w = 0; w < nt / 2; w++)
						{
							int y[3];
							y[0] = 0 + x[0];
							y[1] = 0 + x[1];
							y[2] = 0 + x[2];
							loops_st[w][l].push_back(lat.CalculateWilsonLoop(y, 1, 2, l + 1, w + 1));
						}
					}
				}
			}
		}
		std::cout << "loops calculated\n";
		double val = 0.0;
		for (int w = 0; w < (ns - 1) / 2; w++)
		{
			for (int l = 0; l <= w; l++)
			{
				double num_elements = ((double)loops_ss[l][w].size());
				// calculate the mean
				for (int i = 0; i < num_elements; i++)
				{
					val = ((double)loops_ss[l][w][i]);
					if (val > 1.0)
					{
						val = 1.0;
					}
					else if (val < -1.0)
					{
						val = -1.0;
					}
					means_ss[l][w] += val;
				}
				means_ss[l][w] /= num_elements;
				// calculate the uncertainty
				double variance = 0.0;
				for (int i = 0; i < num_elements; i++)
				{
					val = ((double)loops_ss[l][w][i]);
					if (val > 1.0)
					{
						val = 1.0;
					}
					else if (val < -1.0)
					{
						val = -1.0;
					}
					variance += (val - means_ss[l][w]) * (val - means_ss[l][w]) / num_elements;
				}
				sigmas_ss[l][w] = sqrt(variance / num_elements);
			}
		}
		for (int w = 0; w < (ns - 1) / 2; w++)
		{
			for (int l = 0; l < nt / 2; l++)
			{
				double num_elements = ((double)loops_st[l][w].size());
				// calculate the mean
				for (int i = 0; i < num_elements; i++)
				{
					val = ((double)loops_st[l][w][i]);
					if (val > 1.0)
					{
						val = 1.0;
					}
					else if (val < -1.0)
					{
						val = -1.0;
					}
					means_st[l][w] += val;
				}
				means_st[l][w] /= num_elements;
				// calculate the uncertainty
				double variance = 0.0;
				for (int i = 0; i < num_elements; i++)
				{
					val = ((double)loops_st[l][w][i]);
					if (val > 1.0)
					{
						val = 1.0;
					}
					else if (val < -1.0)
					{
						val = -1.0;
					}
					variance += (val - means_st[l][w]) * (val - means_st[l][w]) / num_elements;
				}
				sigmas_st[l][w] = sqrt(variance / num_elements);
			}
		}
		std::cout << "saving data\n";
		std::ofstream outfile;
		int counter = 0;
		for (int w = 0; w < (ns - 1) / 2; w++)
		{
			for (int l = 0; l <= w; l++)
			{
				outfile.open(wilson_ss_filenames[l][w], std::ios::app);
				outfile << means_ss[l][w] << "," << sigmas_ss[l][w] << "\n";
				outfile.close();
			}
			for (int l = 0; l < nt / 2; l++)
			{
				//counter += 1;
				outfile.open(wilson_st_filenames[w][l], std::ios::app);
				outfile << means_st[l][w] << "," << sigmas_st[l][w] << "\n";
				outfile.close();
			}
		}
		std::cout << "data saved\n";
	}
	if (1 == 1)
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
						std::cout << t1 << t2 << "\n";
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
	}
}
