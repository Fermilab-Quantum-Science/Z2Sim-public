#include "Lattice.h"
#include "WilsonLoopCalc.h"

Lattice::Lattice(int spacial_dim, int temporal_dim, double coup, bool randomstart)
{
	// store the parameters of the lattice
	ns = spacial_dim;
	nt = temporal_dim;
	dim = ns * ns * nt * 3;
	coupling = coup;
	std::cout << coupling << "\n";
	// allocate data for the index array and lattice
	indexarray = new unsigned int[dim];
	link = new short int[dim];
	//wilson_ss = new double[ns * (ns - 1) / 2];
	//wilson_st = new double[(ns - 1) * nt / 2];
	//norms_ss = new double[ns * (ns - 1) / 2];
	//norms_st = new double[(ns - 1) * nt / 2];

	if (randomstart)
	{
		for (int i = 0; i < dim; i++)
		{
			indexarray[i] = i;
			link[i] = 2 * (rand() % 2) - 1;
		}
	}
	else
	{
		for (int i = 0; i < dim; i++)
		{
			indexarray[i] = i;
			link[i] = 1;
		}
	}
}

Lattice::~Lattice()
{
	delete[] link;
	delete[] indexarray;
	delete[] norms_ss;
	delete[] norms_st;
	delete[] wilson_ss;
	delete[] wilson_st;
}

void Lattice::MoveUp(int x[], int d)
{
	x[d] += 1;
	if (d == 2) {
		if (x[d] >= nt) x[d] -= nt;
	}
	return;
}

void Lattice::MoveDown(int x[], int d)
{
	x[d] -= 1;
	if (d == 2) {
		if (x[d] < 0) x[d] += nt;
	}
	return;
}

int Lattice::GetIndex(int x[], int d)
{
	return x[0] + x[1] * ns + x[2] * ns * ns + d * ns * ns * nt;
}

void Lattice::Shuffle() {
	/* shuffle the elements of the index array*/
	if (dim > 1) {
		size_t i;
		for (i = 0; i < dim - 1; i++) {
			size_t j = i + rand() / (RAND_MAX / (dim - i) + 1);
			int t = indexarray[j];
			indexarray[j] = indexarray[i];
			indexarray[i] = t;
		}
	}
}

void Lattice::Update(double betat)
{
	int x[3], d, dperp, staple, staplespacesum, stapletimesum;
	int index;
	double bplus, bminus;
	double beta_s = betat / coupling;
	double beta_t = betat * coupling;
	Shuffle();
	for (int i = 0; i < dim; i++)
	{
		x[0] = indexarray[i] % ns;
		x[1] = (indexarray[i] / ns) % ns;
		x[2] = (indexarray[i] / ns / ns) % nt;
		d = (indexarray[i] / ns / ns / nt) % 3;
		// check if the current site is not an allowable update link
		// if in x direction and left most end of x direction 
		if (d == 0 && x[0] == ns - 1)
		{
			continue;
		}
		// if in y direction and top most end of y direction
		else if (d == 1 && x[1] == ns - 1)
		{
			continue;
		}
		// if at the farthest corner and pointing in spacial direction
		else if ((x[0] == ns - 1 && x[1] == ns - 1) && (d == 0 || d == 1))
		{
			continue;
		}
		// otherwise we are good to go
		else
		{
			staple = 0;
			stapletimesum = 0;
			staplespacesum = 0;
			if (d == 0)
			{
				if (x[1] > 0) {
					MoveDown(x, 1);
					index = GetIndex(x, 1);
					staple = 0 + link[index];
					index = GetIndex(x, 0);
					staple *= link[index];
					MoveUp(x, 0);
					index = GetIndex(x, 1);
					staple *= link[index];
					MoveUp(x, 1);
					MoveDown(x, 0);
					staplespacesum += staple;
				}
				if (x[1] < ns - 1)
				{
					index = GetIndex(x, 1);
					staple = 0 + link[index];
					MoveUp(x, 1);
					index = GetIndex(x, 0);
					staple *= link[index];
					MoveUp(x, 0);
					MoveDown(x, 1);
					index = GetIndex(x, 1);
					MoveDown(x, 0);
					staple *= link[index];
					staplespacesum += staple;
				}
			}
			if (d == 1)
			{
				if (x[0] > 0) {
					MoveDown(x, 0);
					index = GetIndex(x, 0);
					staple = 0 + link[index];
					index = GetIndex(x, 1);
					staple *= link[index];
					MoveUp(x, 1);
					index = GetIndex(x, 0);
					staple *= link[index];
					MoveUp(x, 0);
					MoveDown(x, 1);
					staplespacesum += staple;
				}
				if (x[0] < ns - 1)
				{
					index = GetIndex(x, 0);
					staple = 0 + link[index];
					MoveUp(x, 0);
					index = GetIndex(x, 1);
					staple *= link[index];
					MoveUp(x, 1);
					MoveDown(x, 0);
					index = GetIndex(x, 0);
					staple *= link[index];
					MoveDown(x, 1);
					staplespacesum += staple;
				}
			}
			if (d == 2)
			{
				for (dperp = 0; dperp < 2; dperp++)
				{
					if (x[dperp] < ns - 1)
					{
						staple = 1 * link[GetIndex(x, dperp)];
						MoveUp(x, dperp);
						staple *= link[GetIndex(x, 2)];
						MoveUp(x, 2);
						MoveDown(x, dperp);
						staple *= link[GetIndex(x, dperp)];
						stapletimesum += staple;
						MoveDown(x, 2);
					}
					if (x[dperp] > 0)
					{
						MoveDown(x, dperp);
						staple = link[GetIndex(x, dperp)] * link[GetIndex(x, 2)];
						MoveUp(x, 2);
						staple *= link[GetIndex(x, dperp)];
						MoveUp(x, dperp);
						MoveDown(x, 2);
						stapletimesum += staple;
					}
				}

			}
			else if (d != 2)
			{
				staple = 1;
				staple *= link[GetIndex(x, 2)];
				MoveUp(x, 2);
				staple *= link[GetIndex(x, d)];
				MoveUp(x, d);
				MoveDown(x, 2);
				staple *= link[GetIndex(x, 2)];
				stapletimesum += staple;
				MoveDown(x, 2);
				staple = 1 * link[GetIndex(x, 2)];
				MoveDown(x, d);
				staple *= link[GetIndex(x, d)];
				staple *= link[GetIndex(x, 2)];
				MoveUp(x, 2);
				stapletimesum += staple;
			}
			bplus = exp(beta_s * ((double)staplespacesum) + beta_t * ((double)stapletimesum));
			bminus = 1.0 / bplus;
			bplus = bplus / (bplus + bminus);
			if (((double)rand()) / ((double)RAND_MAX) < bplus) {
				link[GetIndex(x, d)] = 1;
			}
			else {
				link[GetIndex(x, d)] = -1;
			}
		}
	}
}

double Lattice::CalcPlaquetteAverage()
{
	int x[3], d;
	int ind[4];
	double retval = 0;
	int plaq = 0;
	double norm = ((double)(ns - 1) * (ns - 1) * nt);
	for (x[0] = 0; x[0] < ns - 1; x[0]++)
	{
		for (x[1] = 0; x[1] < ns - 1; x[1]++)
		{
			for (x[2] = 0; x[2] < nt; x[2]++)
			{
				ind[0] = x[0] + x[1] * ns + x[2] * ns * ns;
				ind[1] = x[0] + x[1] * ns + x[2] * ns * ns + ns * ns * nt;
				ind[2] = x[0] + (x[1] + 1) * ns + x[2] * ns * ns;
				ind[3] = x[0] + 1 + x[1] * ns + x[2] * ns * ns + ns * ns * nt;
				plaq = link[ind[0]] * link[ind[1]] * link[ind[2]] * link[ind[3]];
				retval += ((double)plaq) / norm;
			}
		}
	}
	return retval;
}

void Lattice::SaveConfiguration(int configuration_number, double betat)
{
	std::string filename = "configuration";
	std::string str1 = std::to_string(betat);
	std::string str2 = std::to_string(coupling);
	str1 = str1.replace(str1.find('.'), 1, "_");
	str2 = str2.replace(str2.find('.'), 1, "_");
	filename += std::to_string(configuration_number);
	filename += "betatau=" + str1;
	filename += "coupling=" + str2;
	filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
	filename += ".csv";
	filename = "../configurations/" + filename;
	std::ofstream outfile;
	outfile.open(filename, std::ios::out | std::ios::trunc);
	for (int i = 0; i < dim; i++)
	{
		outfile << (link[i] + 1) / 2 << "\n";
	}
	outfile.close();

}

void Lattice::LoadConfiguration(int configuration_number, double betat)
{
	std::string filename = "configuration";
	std::string str1 = std::to_string(betat);
	std::string str2 = std::to_string(coupling);
	str1 = str1.replace(str1.find('.'), 1, "_");
	str2 = str2.replace(str2.find('.'), 1, "_");
	filename += std::to_string(configuration_number);
	filename += "betatau=" + str1;
	filename += "coupling=" + str2;
	filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
	filename += ".csv";
	filename = "../configurations/" + filename;
	std::ifstream infile;
	infile.open(filename, std::ifstream::in);
	std::cout << filename << "\n";
	std::string line;
	int counter = 0;
	if (infile.is_open())
	{
		while (std::getline(infile, line))
		{
			link[counter] = 2 * std::stoi(line, nullptr) - 1;
			counter += 1;
		}
	}
}


void Lattice::ResetWilsonLoops()
{
	for (int x = 0; x < ns * (ns + 1) / 2; x++)
	{
		wilson_ss[x] = 0.0;
	}
	for (int x = 0; x < ns * nt / 2; x++)
	{
		wilson_st[x] = 0.0;
	}
}

void Lattice::CalculateNorms()
{
	for (int x = 0; x < ((int)ns - 1) / 2; x++)
	{
		for (int y = 0; y < ((int)ns - 1) / 2; y++)
		{
			if (x > y)
			{
				norms_ss[x + (ns - 1) * y / 2] += 1.0;
			}
			else
			{
				norms_ss[y + (ns - 1) * x / 2] += 1.0;
			}
		}
	}
	for (int x = 0; x < ((int)ns - 1) / 2; x++)
	{
		for (int t = 0; t < nt / 2; t++)
		{
			norms_st[x * (ns - 1) / 2 + t] += 2.0;
		}
	}
}

std::vector<std::vector<double>> Lattice::CalculateOperators(int t)
{
	std::vector<std::vector<std::vector<int>>> operators_at_site;
	std::vector<std::vector<double>> operators_from_edge;
	operators_from_edge.resize(4);
	operators_at_site.resize(4);
	// allocate memory for the array
	for (int i = 0; i < 4; i++)
	{
		operators_from_edge[i].resize((ns - 1) / 2);
		operators_at_site[i].resize((ns - 1));
		for (int j = 0; j < ns - 1; j++)
		{
			operators_at_site[i][j].resize((ns - 1));
		}
	}
	// iterate through all the sites
	for (int x = 0; x < ns - 1; x++)
	{
		for (int y = 0; y < ns - 1; y++)
		{
			int arr[4];
			int v[3];
			for (int i = 0; i < 4; i++)
			{
				v[0] = (x + i % 2) % ns;
				v[1] = (y + i / 2) % ns;
				v[2] = t;
				arr[i] = CalculateWilsonLoop(v, 0, 1, 1, 1);
			}
			operators_at_site[0][x][y] = 0 + arr[0];
			if (x < ns - 2 && y < ns - 2)
			{
				int r01 = arr[0] * arr[1];
				int r02 = arr[0] * arr[2];
				int r23 = arr[2] * arr[3];
				int r13 = arr[1] * arr[3];
				operators_at_site[1][x][y] = 0 + r01 + r02 + r23 + r13;
				operators_at_site[2][x][y] = 0 + arr[0] * (r13 + r23) + arr[2] * (r13 + r01);
				operators_at_site[3][x][y] = r13 * r02;
			}
		}
	}
	//for (int x = 0; x < ns - 1; x++)
	//{
	//	for (int y = 0; y < ns - 1; y++)
	//	{

	//		int arr[9];
	//		int v[3];
	//		for (int i = 0; i < 9; i++)
	//		{
	//			v[0] = (x + i % 3) % ns;
	//			v[1] = (y + i / 3) % ns;
	//			v[2] = t;
	//			arr[i] = CalculateWilsonLoop(v, 0, 1, 1, 1);
	//		}
	//		//=========================
	//		// 1 x 1 plaquette operator
	//		//=========================
	//		operators_at_site[0][x][y] = arr[0];
	//		if (x < ns - 2 && y < ns - 2)
	//		{
	//			//=========================
	//			// 2 x 1 plaquette operator
	//			//=========================
	//			int r01 = arr[0] * arr[1];
	//			int r03 = arr[0] * arr[3];
	//			int r14 = arr[1] * arr[4];
	//			int r34 = arr[3] * arr[4];
	//			operators_at_site[1][x][y] = r01 + r03 + r14 + r34;
	//			std::cout << operators_at_site[1][x][y];
	//			//=========================
	//			// x0   0x   xx   xx
	//			// xx + xx + x0 + 0x
	//			//=========================
	//			operators_at_site[2][x][y] = arr[1] * (r03 + r34) + arr[3] * (r01 + r14);
	//			//=========================
	//			// 2x2 plaquette operator
	//			//=========================
	//			operators_at_site[3][x][y] = r01 * r34;
	//			//if (x < ns - 3 && y < ns - 3)
	//			//{
	//			//	//=========================
	//			//	// 3 x 1 rectangle
	//			//	//=========================
	//			//	int r147 = arr[1] * arr[4] * arr[7];
	//			//	int r345 = arr[3] * arr[4] * arr[5];
	//			//	operators_at_site[4][x][y] += r147 + r345;
	//			//	//=========================
	//			//	// 2x3 L operators 
	//			//	//=========================
	//			//	int corners = (arr[2] + arr[8] + arr[0] + arr[6]);
	//			//	operators_at_site[5][x][y] += corners * (r345 + r147);
	//			//	//=========================
	//			//	//2x3 T operators
	//			//	//=========================
	//			//	int t = r345 * (arr[1] + arr[7]);
	//			//	operators_at_site[6][x][y] += t + r147 * (arr[5] + arr[3]);
	//			//	//=========================
	//			//	// 2x3 P operators
	//			//	//=========================
	//			//	operators_at_site[7][x][y] = r345 * ((arr[6] + arr[8]) * arr[7] + (arr[0] + arr[2]) * arr[1]);
	//			//	operators_at_site[7][x][y] += r147 * ((arr[6] + arr[0]) * arr[3] + (arr[2] + arr[8]) * arr[5]);
	//			//	//=========================
	//			//	// 2x3 U operator
	//			//	//=========================
	//			//	int U_up = r147 * arr[2] * arr[8];
	//			//	int U_down = r147 * arr[0] * arr[6];
	//			//	int U_right = r345 * arr[6] * arr[8];
	//			//	int U_left = r345 * arr[0] * arr[2];
	//			//	operators_at_site[8][x][y] = U_up + U_down + U_left + U_right;
	//			//	//=========================
	//			//	// 2x3 rectangles
	//			//	//=========================
	//			//	int twobythree = r345 * arr[0] * arr[1] * arr[2];
	//			//	int threebytwo = r147 * arr[2] * arr[5] * arr[8];
	//			//	operators_at_site[9][x][y] = twobythree + threebytwo;
	//			//	//=========================
	//			//	// 3x3 T
	//			//	//=========================
	//			//	operators_at_site[10][x][y] += r345 * arr[1] * arr[7];
	//			//	//=========================
	//			//	// 3x3 Ls
	//			//	//=========================
	//			//	int set1 = arr[0] * arr[8];
	//			//	int set2 = arr[2] * arr[6];
	//			//	int LDL = set2 * arr[1] * arr[0] * arr[3];
	//			//	int LUL = set1 * arr[1] * arr[2] * arr[5];
	//			//	int LUR = set2 * arr[5] * arr[8] * arr[7];
	//			//	int LDR = set1 * arr[3] * arr[7] * arr[6];
	//			//	int LSET = LUL + LDL + LUR + LDR;
	//			//	operators_at_site[11][x][y] += LSET;
	//			//	//=========================
	//			//	// 3x3 J's
	//			//	//=========================
	//			//	operators_at_site[12][x][y] = LDL * (arr[7] + arr[5]);
	//			//	operators_at_site[12][x][y] += LUL * (arr[3] + arr[7]);
	//			//	operators_at_site[12][x][y] += LUR * (arr[1] + arr[3]);
	//			//	operators_at_site[12][x][y] += LDR * (arr[1] + arr[5]);
	//			//	//=========================
	//			//	// 3x3 step
	//			//	//==========================
	//			//	operators_at_site[13][x][y] += LSET * arr[4];
	//			//	//=========================
	//			//	// 3x3 edge 
	//			//	// x00
	//			//	// xxx
	//			//	// xxx
	//			//	//=========================
	//			//	int P1 = twobythree * (arr[8] + arr[6]);
	//			//	int P2 = threebytwo * (arr[0] + arr[6]);
	//			//	int P4 = r147 * arr[0] * arr[3] * arr[6] * (arr[8] + arr[2]);
	//			//	int P3 = r345 * arr[6] * arr[7] * arr[8] * (arr[0] + arr[2]);
	//			//	operators_at_site[14][x][y] += P1 + P2 + P4 + P3;
	//			//	//=========================
	//			//	// 3x3 U's
	//			//	//=========================
	//			//	int fourcorn = arr[0] * arr[2] * arr[6] * arr[8];
	//			//	int U3x3 = fourcorn * (arr[1] * arr[3] * (arr[7] + arr[5]) + (arr[1] + arr[3]) * arr[7] * arr[5]);
	//			//	operators_at_site[15][x][y] += U3x3;
	//			//	//=========================
	//			//	// 3x3 Chonky U's
	//			//	//=========================
	//			//	operators_at_site[16][x][y] += U3x3 * arr[4];
	//			//	//=========================
	//			//	// 3x3 T's
	//			//	//=========================
	//			//	int T1 = r345 * arr[0] * arr[6];
	//			//	int T2 = r345 * arr[2] * arr[8];
	//			//	int T3 = r147 * arr[0] * arr[2];
	//			//	int T4 = r147 * arr[6] * arr[8];
	//			//	operators_at_site[17][x][y] = T1 + T2 + T3 + T4;
	//			//	//=========================
	//			//	// 3x3 lowercase h
	//			//	//=========================
	//			//	operators_at_site[18][x][y] = T1 * (arr[2] + arr[8]) + T2 * (arr[0] + arr[6]);
	//			//	operators_at_site[18][x][y] += T3 * (arr[8] + arr[6]) + T4 * (arr[0] + arr[2]);
	//			//	//=========================
	//			//	// 3x3 blob
	//			//	//=========================
	//			//	operators_at_site[19][x][y] = T1 * (arr[2] * arr[7] + arr[1] * arr[8]) + T2 * (arr[0] * arr[7] + arr[1] * arr[6]);
	//			//	operators_at_site[19][x][y] += T3 * (arr[8] * arr[3] + arr[6] * arr[5]) + T4 * (arr[0] * arr[5] + arr[2] * arr[3]);
	//			//	//=========================
	//			//	// 3x3 S
	//			//	//=========================
	//			//	operators_at_site[20][x][y] = (r147 + r345) * (arr[0] * arr[8] + arr[2] * arr[6]);
	//			//	//=========================
	//			//	// 3x3 blob 2
	//			//	//=========================
	//			//	operators_at_site[21][x][y] = (r147 * (arr[3] + arr[5]) + r345 * (arr[1] + arr[7])) * (arr[0] * arr[8] + arr[2] * arr[6]);
	//			//	//=========================
	//			//	// 3x3 Tblob
	//			//	//=========================
	//			//	operators_at_site[22][x][y] = (T2 + T1) * (arr[1] + arr[7]) + (T3 + T4) * (arr[3] + arr[5]);
	//			//	//=========================
	//			//	// 3x3 Chonky T
	//			//	//=========================
	//			//	operators_at_site[23][x][y] = r345 * arr[1] * arr[7] * (arr[0] * arr[2] + arr[8] * arr[6]);
	//			//	operators_at_site[23][x][y] += r147 * arr[3] * arr[5] * (arr[0] * arr[6] + arr[2] * arr[8]);
	//			//	//=========================
	//			//	// 3x3 Capital H
	//			//	//=========================
	//			//	operators_at_site[24][x][y] = (r147 + r345) * fourcorn;
	//			//	//=========================
	//			//	// 3x3 Slopy
	//			//	//=========================
	//			//	operators_at_site[25][x][y] = r147 * ((arr[0] + arr[6]) * arr[5] + (arr[8] + arr[2]) * arr[3]);
	//			//	operators_at_site[25][x][y] += r345 * (arr[7] * (arr[0] + arr[2]) + arr[1] * (arr[6] + arr[8]));
	//			//	//=========================
	//			//	// 3x3 square
	//			//	//=========================
	//			//	int threebythree = threebytwo * arr[0] * arr[3] * arr[5];
	//			//	operators_at_site[27][x][y] = threebythree;
	//			//	//=========================
	//			//	// trunc 3x3
	//			//	//=========================
	//			//	operators_at_site[26][x][y] = threebythree * corners;
	//			// }
	//		}
	//	}
	//}
	for (int d = 0; d < (ns - 1) / 2; d++)
	{
		for (int x = d; x < ns - 1 - d; x++)
		{
			for (int y = d; y < ns - 1 - d; y++)
			{
				operators_from_edge[0][d] += ((double)operators_at_site[0][x][y]);
				if ((ns - 1) / 2 - d > 0)
				{
					for (int j = 1; j < 4; j++)
					{
						operators_from_edge[j][d] += ((double)operators_at_site[j][x][y]);
					}
				}
				/*if ((ns - 1) / 2 - d> 1)
				{
					for (int j = 4; j < 28; j++)
					{
						operators_from_edge[j][d] += ((double)operators_at_site[j][x][y]);
					}
				}*/
			}
		}
	}
	for (int d = 0; d < (ns - 1) / 2; d++)
	{
		double norm = ((double)(ns - 1 - d) * (ns - 1 - d));
		for (int i = 0; i < 4; i++)
		{
			operators_from_edge[i][d] /= sqrt(norm);
		}
	}
	return operators_from_edge;
}

void Lattice::CalculateCrossCorrelatorsOnTheFly(int ns, int nt, double coupling, double betatau, std::vector<std::vector<std::string>> operator_filenames, std::vector<std::vector<std::vector<std::vector<std::string>>>> correlatorfilenames)
{
	int x[3];
	int opernum = 4;
	std::cout << "Calculating operators\n";
	std::vector<std::vector<std::vector<double> > > operators;
	std::vector<std::vector<std::vector<std::vector<double> > > > correlators;
	operators.resize(opernum);
	correlators.resize(opernum);
	for (int i = 0; i < opernum; i++)
	{
		operators[i].resize(nt);
		for (int j = 0; j < nt; j++)
		{
			operators[i][j].resize((ns - 1) / 2);
		}
		correlators[i].resize(opernum);
		for (int j = 0; j < opernum; j++)
		{
			correlators[i][j].resize(nt / 2);
			for (int t = 0; t < nt / 2; t++)
			{
				correlators[i][j][t].resize((ns - 1) / 2);
			}
		}
	}
	for (int t = 0; t < nt; t++)
	{
		std::vector<std::vector<double> > operators_1 = CalculateOperators(t);
		for (int op1 = 0; op1 < opernum; op1++)
		{
			for (int d = 0; d < (ns - 1) / 2; d++)
			{
				operators[op1][t][d] += operators_1[op1][d];
			}
		}
	}
	for (int d = 0; d < (ns - 1) / 2; d++)
	{
		for (int t = 0; t < nt; t++)
		{
			for (int t0 = 0; t0 < nt / 2; t0++)
			{
				for (int op1 = 0; op1 < opernum; op1++)
				{
					for (int op2 = 0; op2 < opernum; op2++)
					{
						correlators[op1][op2][t0][d] += operators[op1][t][d] * operators[op2][(t + t0) % nt][d] / ((double)nt);
					}
				}
			}
		}
	}
	std::ofstream outfile;
	for (int d = 0; d < (ns - 1) / 2; d++)
	{
		for (int op1 = 0; op1 < opernum; op1++)
		{
			float opval = 0.0;
			for (int t = 0; t < nt; t++)
			{
				opval += operators[op1][t][d] / ((float)nt);
			}
			outfile.open(operator_filenames[op1][d], std::ios::app);
			outfile << opval << "\n";
			outfile.close();
			for (int op2 = 0; op2 < opernum; op2++)
			{
				for (int dt = 0; dt < nt / 2; dt++)
				{
					outfile.open(correlatorfilenames[op1][op2][dt][d], std::ios::app);
					outfile << correlators[op1][op2][dt][d] << ", 0.0\n";
					outfile.close();
				}
			}
		}
	}
}

void Lattice::CalculateWilsonLoopsLocalandCrossCorrelatorsOnTheFly(int ns, int nt, double coupling, double betatau, int montecarlosweeps, int configurations)
{
	std::string filename;
	// [x-length of loop][y-length of loop][distance from edge][configuration]
	std::vector<std::vector<std::vector<std::vector<float> > > > ss_loops;
	// [x-length of loop][t-length of loop][distance from edge][configuration]
	std::vector<std::vector<std::vector<std::vector<float> > > > st_loops;
	// [operator][t-length][distance from edge][configuration]
	std::vector<std::vector<std::vector<std::vector<float> > > > opers;
	// [operator][distance from edge][configuration]
	std::vector<std::vector<std::vector<float> > > average;
	// holder vector
	// [operator][time][distance from edge]
	std::vector<std::vector<std::vector<float> > > holder;
	holder.resize(4);
	ss_loops.resize((ns - 1) / 2);
	st_loops.resize((ns - 1) / 2);
	opers.resize(4);
	average.resize(4);

	std::vector<std::vector<std::vector<float> > > st_loops_norm;
	std::vector<std::vector<std::vector< int > > > ss_loops_norm;
	ss_loops.resize((ns - 1) / 2);
	st_loops.resize((ns - 1) / 2);
	st_loops_norm.resize((ns - 1) / 2);
	ss_loops_norm.resize((ns - 1) / 2);
	for (int x = 0; x < (ns - 1) / 2; x++)
	{
		st_loops_norm[x].resize(nt / 2);
		ss_loops_norm[x].resize((ns - 1) / 2);
		for (int y = 0; y < (ns - 1) / 2; y++)
		{
			ss_loops_norm[x][y].resize((ns - 1) / 2);
		}
		for (int y = 0; y < nt / 2; y++)
		{
			st_loops_norm[x][y].resize((ns - 1) / 2);
		}
	}
	int v[3];
	bool con1, con2, con3, con4;
	for (int l = 1; l <= (ns - 1) / 2; l++)
	{
		for (int w = 1; w <= (ns - 1) / 2; w++)
		{
			for (v[2] = 0; v[2] < nt; v[2]++)
			{
				for (v[0] = 0; v[0] < ns - l; v[0]++)
				{
					for (v[1] = 0; v[1] < ns - w; v[1]++)
					{
						for (int d = 0; d < (ns - 1) / 2; d++)
						{
							con1 = (d <= v[0]);
							con2 = (d <= v[1]);
							con3 = (v[0] + l < ns - d);
							con4 = (v[1] + w < ns - d);
							if (con1 && con2 && con3 && con4)
							{
								ss_loops_norm[l - 1][w - 1][d] += 1;
							}
						}
					}
				}
			}

		}
	}
	for (int l = 1; l <= (ns - 1) / 2; l++)
	{
		for (int w = 1; w <= nt / 2; w++)
		{
			for (v[2] = 0; v[2] < nt; v[2]++)
			{
				for (v[0] = 0; v[0] < ns - l; v[0]++)
				{
					for (v[1] = 0; v[1] < ns - 1; v[1]++)
					{
						for (int d = 0; d < (ns - 1) / 2; d++)
						{
							//std::cout << v[0] << "," << v[1] << "," << v[2] << "," << d << "\n";
							con1 = (d <= v[0]);
							con2 = (d <= v[1]);
							con3 = (v[0] + l < ns - d);
							con4 = (v[1] < ns - d);
							if (con1 && con2 && con3 && con4)
							{
								st_loops_norm[l - 1][w - 1][d] += 1.0;
							}
						}
					}
				}
				for (v[0] = 0; v[0] < ns - 1; v[0]++)
				{
					for (v[1] = 0; v[1] < ns - l; v[1]++)
					{
						for (int d = 0; d < (ns - 1) / 2; d++)
						{
							con1 = (d <= v[0]);
							con2 = (d <= v[1]);
							con3 = (v[0] < ns - d);
							con4 = (v[1] + l < ns - d);
							if (con1 && con2 && con3 && con4)
							{
								st_loops_norm[l - 1][w - 1][d] += 1;
							}
						}
					}
				}
			}
		}
	}
	int arrsize = 1000;
	for (int op = 0; op < 4; op++)
	{
		holder[op].resize(nt);
		for (int t = 0; t < nt; t++)
		{
			holder[op][t].resize((ns - 1) / 2);
		}

		opers[op].resize(nt / 2);
		average[op].resize((ns - 1) / 2);
		for (int d = 0; d < (ns - 1) / 2; d++)
		{
			average[op][d].resize(arrsize);
		}
		for (int t = 0; t < nt / 2; t++)
		{
			opers[op][t].resize((ns - 1) / 2);
			for (int d = 0; d < (ns - 1) / 2; d++)
			{
				opers[op][t][d].resize(arrsize);
			}
		}
	}
	for (int x = 0; x < (ns - 1) / 2; x++)
	{
		ss_loops[x].resize((ns - 1) / 2);
		st_loops[x].resize(nt / 2);
		for (int y = 0; y < (ns - 1) / 2; y++)
		{
			ss_loops[x][y].resize((ns - 1) / 2);
			for (int d = 0; d < (ns - 1) / 2; d++)
			{
				ss_loops[x][y][d].resize(arrsize);
			}
		}
		for (int t = 0; t < nt / 2; t++)
		{
			st_loops[x][t].resize((ns - 1) / 2);
			for (int d = 0; d < (ns - 1) / 2; d++)
			{
				st_loops[x][t][d].resize(arrsize);
			}
		}
	}

	for (int config = 0; config < configurations; config++)
	{
		for (int i = 0; i < montecarlosweeps; i++)
		{
			Update(betatau);
		}
		std::cout << "configuration" << config << "generated\n";
		short wloop = 0;
		bool con1, con2, con3, con4;
		short norm = 0;
		for (int t = 0; t < nt; t++)
		{
			std::vector<std::vector<double> > hold = CalculateOperators(t);
			for (int op = 0; op < 4; op++)
			{
				for (int d = 0; d < (ns - 1) / 2; d++)
				{
					holder[op][t][d] = 0.0 + hold[op][d];
					average[op][d][config % arrsize] += hold[op][d] / ((double)nt);
				}
			}
		}
		std::cout << "correlators calculated\n";
		for (int t = 0; t < nt; t++)
		{
			for (int dt = 0; dt < nt / 2; dt++)
			{
				for (int op = 0; op < 4; op++)
				{
					for (int d = 0; d < (ns - 1) / 2; d++)
					{
						opers[op][dt][d][config % arrsize] += ((float) holder[op][t][d] * holder[op][(t + dt) % nt][d]) / ((float)nt);
					}
				}
			}
		}
		std::cout << "correlators stored\n";
		/*
		for (int l = 1; l <= (ns - 1) / 2; l++)
		{
			for (int w = 1; w <= (ns - 1) / 2; w++)
			{
				for (v[2] = 0; v[2] < nt; v[2]++)
				{
					for (v[0] = 0; v[0] < ns - l; v[0]++)
					{
						for (v[1] = 0; v[1] < ns - w; v[1]++)
						{
							for (int d = 0; d < (ns - 1) / 2; d++)
							{
								//std::cout << v[0] << "," << v[1] << "," << v[2] << "," << d << "\n";
								wloop = ((float)CalculateWilsonLoop(v, 0, 1, l, w));
								con1 = (d <= v[0]);
								con2 = (d <= v[1]);
								con3 = (v[0] + l < ns - d);
								con4 = (v[1] + w < ns - d);
								if (con1 && con2 && con3 && con4)
								{
									ss_loops[l - 1][w - 1][d][config % arrsize] += wloop;
								}
							}
						}
					}
				}

			}
		}
		for (int l = 1; l <= (ns - 1) / 2; l++)
		{
			for (int w = 1; w <= nt / 2; w++)
			{
				for (v[2] = 0; v[2] < nt; v[2]++)
				{
					for (v[0] = 0; v[0] < ns - l; v[0]++)
					{
						for (v[1] = 0; v[1] < ns - w; v[1]++)
						{
							wloop = ((float)CalculateWilsonLoop(v, 0, 2, l, w));
							for (int d = 0; d < (ns - 1) / 2; d++)
							{
								//std::cout << v[0] << "," << v[1] << "," << v[2] << "," << d << "\n";
								con1 = (d <= v[0]);
								con2 = (d <= v[1]);
								con3 = (v[0] + l < ns - d);
								con4 = (v[1] < ns - d);
								if (con1 && con2 && con3 && con4)
								{
									st_loops[l - 1][w - 1][d][config % arrsize] += wloop;
								}
							}
							wloop = ((float)CalculateWilsonLoop(v, 1, 2, l, w));
							for (int d = 0; d < (ns - 1) / 2; d++)
							{
								con1 = (d <= v[0]);
								con2 = (d <= v[1]);
								con3 = (v[0] < ns - d);
								con4 = (v[1] + l < ns - d);
								if (con1 && con2 && con3 && con4)
								{
									st_loops[l - 1][w - 1][d][config % arrsize] += wloop;
								}
							}
						}
					}
				}
			}
		}
		*/
		//short wloop = 0;
		//bool con1, con2, con3, con4;
		//short norm = 0;
		for (int l = 1; l <= (ns - 1) / 2; l++)
		{
			for (int w = 1; w <= (ns - 1) / 2; w++)
			{
				for (v[2] = 0; v[2] < nt; v[2]++)
				{
					for (v[0] = 0; v[0] < ns - l; v[0]++)
					{
						for (v[1] = 0; v[1] < ns - w; v[1]++)
						{
							for (int d = 0; d < (ns - 1) / 2; d++)
							{
								//std::cout << v[0] << "," << v[1] << "," << v[2] << "," << d << "\n";
								wloop = ((float)CalculateWilsonLoop(v, 0, 1, l, w));
								con1 = (d <= v[0]);
								con2 = (d <= v[1]);
								con3 = (v[0] + l < ns - d);
								con4 = (v[1] + w < ns - d);
								if (con1 && con2 && con3 && con4)
								{
									ss_loops[l - 1][w - 1][d][config % arrsize] += wloop;
								}
							}
						}
					}
				}

			}
		}
		for (int l = 1; l <= (ns - 1) / 2; l++)
		{
			for (int w = 1; w <= nt / 2; w++)
			{
				for (v[2] = 0; v[2] < nt; v[2]++)
				{
					for (v[0] = 0; v[0] < ns - l; v[0]++)
					{
						for (v[1] = 0; v[1] < ns - 1; v[1]++)
						{
							wloop = ((float)CalculateWilsonLoop(v, 0, 2, l, w));
							for (int d = 0; d < (ns - 1) / 2; d++)
							{
								//std::cout << v[0] << "," << v[1] << "," << v[2] << "," << d << "\n";
								con1 = (d <= v[0]);
								con2 = (d <= v[1]);
								con3 = (v[0] + l < ns - d);
								con4 = (v[1] < ns - d);
								if (con1 && con2 && con3 && con4)
								{
									st_loops[l - 1][w - 1][d][config % arrsize] += wloop;
								}
							}
						}
					}
					for (v[0] = 0; v[0] < ns - 1; v[0]++)
					{
						for (v[1] = 0; v[1] < ns - l; v[1]++)
						{
							wloop = ((float)CalculateWilsonLoop(v, 1, 2, l, w));
							for (int d = 0; d < (ns - 1) / 2; d++)
							{
								con1 = (d <= v[0]);
								con2 = (d <= v[1]);
								con3 = (v[0] < ns - d);
								con4 = (v[1] + l < ns - d);
								if (con1 && con2 && con3 && con4)
								{
									st_loops[l - 1][w - 1][d][config % arrsize] += wloop;
								}
							}
						}
					}
				}
			}
		}
		std::cout << "wilson loops calculated\n";
		if ((config + 1) % arrsize == 0)
		{

			std::ofstream outfile;
			std::string str1 = std::to_string(betatau);
			std::string str2 = std::to_string(coupling);
			str1 = str1.replace(str1.find('.'), 1, "_");
			str2 = str2.replace(str2.find('.'), 1, "_");

			for (int op = 0; op < 4; op++)
			{
				for (int d = 0; d < (ns - 1) / 2; d++)
				{
					for (int dt = 0; dt < nt / 2; dt++)
					{
						filename = "../correlators/";
						filename += "correlator=" + std::to_string(op);
						filename += "rim=" + std::to_string(d) + "dt=" + std::to_string(dt) + "beta=" + str1;
						filename += "anisotropy=" + str2;
						filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
						filename += ".csv";
						outfile.open(filename, std::ios::app);
						std::cout << "openning file\n";
						for (int c = 0; c < arrsize; c++)
						{
							outfile << opers[op][dt][d][c] << "\n";
							opers[op][dt][d][c] = 0;
						}
						outfile.close();
					}
					filename = "../correlators/";
					filename += "operator=" + std::to_string(op);
					filename += "rim=" + std::to_string(d) + "beta=" + str1;
					filename += "anisotropy=" + str2;
					filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
					filename += ".csv";
					outfile.open(filename, std::ios::app);
					for (int c = 0; c < arrsize; c++)
					{
						outfile << average[op][d][c] << "\n";
						average[op][d][c] = 0.0;
					}
					outfile.close();
				}
			}
			std::cout << "Correlators written too\n";
			for (int l = 0; l < (ns - 1) / 2; l++)
			{
				for (int w = 0; w < (ns - 1) / 2; w++)
				{
					filename = "../wilsonloops/";
					filename += "wilsonloops_ss_l=" + std::to_string(l) + "w=" + std::to_string(w) + "beta=" + str1;
					filename += "anisotropy=" + str2;
					filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
					filename += ".csv";
					outfile.open(filename, std::ios::app);
					for (int c = 0; c < arrsize; c++)
					{
						int d = 0;
						for (d = 0; d < (ns - 1) / 2 - 1; d++)
						{
							outfile << ss_loops[l][w][d][c] / ((float)ss_loops_norm[l][w][d]) << ",";// ((float)std::max(1, (ns - l - 1 - 2 * d)* (ns - w - 1 - 2 * d)* nt)) << ",";
							ss_loops[l][w][d][c] = 0;
						}
						outfile << ss_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)ss_loops_norm[l][w][(ns - 1) / 2 - 1]) << "\n";// ((float)(ns - l - 1 - 2 * (d))* (ns - w - 1 - 2 * (d))* nt) << "\n";
						//outfile << ss_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)std::max(1, (ns - l - 1 - 2 * (d + 1)) * (ns - w - 1 - 2 * (d + 1)) * nt)) << "\n";
						ss_loops[l][w][(ns - 1) / 2 - 1][c] = 0;
					}
					outfile.close();
				}
				for (int w = 0; w < nt / 2; w++)
				{
					filename = "../wilsonloops/";
					filename += "wilsonloops_st_l=" + std::to_string(l) + "w=" + std::to_string(w) + "beta=" + str1;
					filename += "anisotropy=" + str2;
					filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
					filename += ".csv";
					outfile.open(filename, std::ios::app);
					for (int c = 0; c < arrsize; c++)
					{
						int d = 0;
						for (d = 0; d < (ns - 1) / 2 - 1; d++)
						{
							float norm = ((float)(ns - l - 2 * d) * (ns - 2 * d) * 2 * nt);
							outfile << st_loops[l][w][d][c] / ((float)st_loops_norm[l][w][d]) << ",";// std::max(1.0f, norm) << ",";
							st_loops[l][w][d][c] = 0;
						}
						//outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)std::max(1, (ns - l - 2 * (d + 1)) * (ns - 2 * d) * 2 * nt)) << "\n";
						outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)st_loops_norm[l][w][(ns - 1) / 2 - 1]) << ",";// / ((float)(ns - l - 2 * (d)) * (ns - 2 * d) * 2 * nt) << "\n";
						//outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)std::max(1, (ns - l - 2 * (d + 1)) * (ns - 2 * (d + 1)) * 2 * nt)) << "\n";
						st_loops[l][w][(ns - 1) / 2 - 1][c] = 0;
					}
					outfile.close();
				}
			}
		}
		std::cout << "files written to\n";
	}
	std::cout << "Saving\n";
	std::ofstream outfile;
	std::string str1 = std::to_string(betatau);
	std::string str2 = std::to_string(coupling);
	str1 = str1.replace(str1.find('.'), 1, "_");
	str2 = str2.replace(str2.find('.'), 1, "_");
	std::cout << str1 << str2 << "\n";
	for (int op = 0; op < 4; op++)
	{
		for (int d = 0; d < (ns - 1) / 2; d++)
		{
			for (int dt = 0; dt < nt / 2; dt++)
			{
				std::cout << op << " " << d << " " << dt << "\n";
				filename = "../correlators/";
				filename += "correlator=" + std::to_string(op);
				filename += "rim=" + std::to_string(d) + "dt=" + std::to_string(dt) + "beta=" + str1;
				filename += "anisotropy=" + str2;
				filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
				filename += ".csv";
				std::cout << configurations % arrsize << "\n";
				outfile.open(filename, std::ios::app);
				for (int c = 0; c < configurations % arrsize; c++)
				{
					outfile << opers[op][dt][d][c] << "\n";
				}
				outfile.close();
			}
			filename = "../correlators/";
			filename += "operator=" + std::to_string(op);
			filename += "rim=" + std::to_string(d) + "beta=" + str1;
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			filename += ".csv";
			outfile.open(filename, std::ios::app);
			for (int c = 0; c < configurations % arrsize; c++)
			{
				outfile << average[op][d][c] << "\n";
			}
			outfile.close();
		}
	}
	for (int l = 0; l < (ns - 1) / 2; l++)
	{
		for (int w = 0; w < (ns - 1) / 2; w++)
		{
			filename = "../wilsonloops/";
			filename += "wilsonloops_ss_l=" + std::to_string(l) + "w=" + std::to_string(w) + "beta=" + str1;
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			filename += ".csv";
			outfile.open(filename, std::ios::app);
			for (int c = 0; c < configurations % arrsize; c++)
			{
				int d = 0;
				for (d = 0; d < (ns - 1) / 2 - 1; d++)
				{
					outfile << ss_loops[l][w][d][c] / ((float)(ns - l - 1 - 2 * d) * (ns - w - 1 - 2 * d) * nt) << ",";
				}
				outfile << ss_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)(ns - l - 1 - 2 * (d)) * (ns - w - 1 - 2 * (d)) * nt) << "\n";
			}
			outfile.close();
		}
		for (int w = 0; w < nt / 2; w++)
		{
			filename = "../wilsonloops/";
			filename += "wilsonloops_st_l=" + std::to_string(l) + "w=" + std::to_string(w) + "beta=" + str1;
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			filename += ".csv";
			outfile.open(filename, std::ios::app);
			for (int c = 0; c < configurations % arrsize; c++)
			{
				int d = 0;
				for (d = 0; d < (ns - 1) / 2 - 1; d++)
				{
					float norm = ((float)(ns - l - 2 * d) * (ns - 2 * d) * 2 * nt);
					outfile << st_loops[l][w][d][c] / norm << ",";
				}
				outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)(ns - l - 2 * (d)) * (ns - 2 * d) * 2 * nt) << "\n";
			}
			outfile.close();
		}
	}
}

void Lattice::CalculateWilsonLoopsLocalOnTheFly(int ns, int nt, double coupling, double betatau, int montecarlosweeps, int configurations)
{
	std::string filename;
	// [x-length of loop][y-length of loop][distance from edge][configuration]
	std::vector<std::vector<std::vector<std::vector<float> > > > ss_loops;
	// [x-length of loop][t-length of loop][distance from edge][configuration]
	std::vector<std::vector<std::vector<std::vector<float> > > > st_loops;
	std::vector<std::vector<std::vector<float> > > st_loops_norm;
	std::vector<std::vector<std::vector< int > > > ss_loops_norm;
	ss_loops.resize((ns - 1) / 2);
	st_loops.resize((ns - 1) / 2);
	st_loops_norm.resize((ns - 1) / 2);
	ss_loops_norm.resize((ns - 1) / 2);
	for (int x = 0; x < (ns - 1) / 2; x++)
	{
		st_loops_norm[x].resize(nt / 2);
		ss_loops_norm[x].resize((ns - 1) / 2);
		for (int y = 0; y < (ns - 1) / 2; y++)
		{
			ss_loops_norm[x][y].resize((ns - 1) / 2);
		}
		for (int y = 0; y < nt / 2; y++)
		{
			st_loops_norm[x][y].resize((ns - 1) / 2);
		}
	}
	int v[3];
	bool con1, con2, con3, con4;
	for (int l = 1; l <= (ns - 1) / 2; l++)
	{
		for (int w = 1; w <= (ns - 1) / 2; w++)
		{
			for (v[2] = 0; v[2] < nt; v[2]++)
			{
				for (v[0] = 0; v[0] < ns - l; v[0]++)
				{
					for (v[1] = 0; v[1] < ns - w; v[1]++)
					{
						for (int d = 0; d < (ns - 1) / 2; d++)
						{
							con1 = (d <= v[0]);
							con2 = (d <= v[1]);
							con3 = (v[0] + l < ns - d);
							con4 = (v[1] + w < ns - d);
							if (con1 && con2 && con3 && con4)
							{
								ss_loops_norm[l - 1][w - 1][d] += 1;
							}
						}
					}
				}
			}

		}
	}
	for (int l = 1; l <= (ns - 1) / 2; l++)
	{
		for (int w = 1; w <= nt / 2; w++)
		{
			for (v[2] = 0; v[2] < nt; v[2]++)
			{
				for (v[0] = 0; v[0] < ns - l; v[0]++)
				{
					for (v[1] = 0; v[1] < ns - 1; v[1]++)
					{
						for (int d = 0; d < (ns - 1) / 2; d++)
						{
							//std::cout << v[0] << "," << v[1] << "," << v[2] << "," << d << "\n";
							con1 = (d <= v[0]);
							con2 = (d <= v[1]);
							con3 = (v[0] + l < ns - d);
							con4 = (v[1] < ns - d);
							if (con1 && con2 && con3 && con4)
							{
								st_loops_norm[l - 1][w - 1][d]+= 1.0;
							}
						}
					}
				}
				for (v[0] = 0; v[0] < ns - 1; v[0]++)
				{
					for (v[1] = 0; v[1] < ns - l; v[1]++)
					{
						for (int d = 0; d < (ns - 1) / 2; d++)
						{
							con1 = (d <= v[0]);
							con2 = (d <= v[1]);
							con3 = (v[0] < ns - d);
							con4 = (v[1] + l < ns - d);
							if (con1 && con2 && con3 && con4)
							{
								st_loops_norm[l - 1][w - 1][d] += 1;
							}
						}
					}
				}
			}
		}
	}
	int arrsize = 1000;
	for (int x = 0; x < (ns - 1) / 2; x++)
	{
		ss_loops[x].resize((ns - 1) / 2);
		st_loops[x].resize(nt / 2);
		for (int y = 0; y < (ns - 1) / 2; y++)
		{
			ss_loops[x][y].resize((ns - 1) / 2);
			for (int d = 0; d < (ns - 1) / 2; d++)
			{
				ss_loops[x][y][d].resize(arrsize);
			}
		}
		for (int t = 0; t < nt / 2; t++)
		{
			st_loops[x][t].resize((ns - 1) / 2);
			for (int d = 0; d < (ns - 1) / 2; d++)
			{
				st_loops[x][t][d].resize(arrsize);
			}
		}
	}
	for (int config = 0; config < configurations; config++)
	{
		for (int i = 0; i < montecarlosweeps; i++)
		{
			Update(betatau);
		}
		std::cout << "configuration" << config << "generated\n";
		short wloop = 0;
		short norm = 0;
		for (int l = 1; l <= (ns - 1) / 2; l++)
		{
			for (int w = 1; w <= (ns - 1) / 2; w++)
			{
				for (v[2] = 0; v[2] < nt; v[2]++)
				{
					for (v[0] = 0; v[0] < ns - l; v[0]++)
					{
						for (v[1] = 0; v[1] < ns - w; v[1]++)
						{
							for (int d = 0; d < (ns - 1) / 2; d++)
							{
								//std::cout << v[0] << "," << v[1] << "," << v[2] << "," << d << "\n";
								wloop = ((float)CalculateWilsonLoop(v, 0, 1, l, w));
								con1 = (d <= v[0]);
								con2 = (d <= v[1]);
								con3 = (v[0] + l < ns - d);
								con4 = (v[1] + w < ns - d);
								if (con1 && con2 && con3 && con4)
								{
									ss_loops[l - 1][w - 1][d][config % arrsize] += wloop;
								}
							}
						}
					}
				}

			}
		}
		for (int l = 1; l <= (ns - 1) / 2; l++)
		{
			for (int w = 1; w <= nt / 2; w++)
			{
				for (v[2] = 0; v[2] < nt; v[2]++)
				{
					for (v[0] = 0; v[0] < ns - l; v[0]++)
					{
						for (v[1] = 0; v[1] < ns - 1; v[1]++)
						{
							wloop = ((float)CalculateWilsonLoop(v, 0, 2, l, w));
							for (int d = 0; d < (ns - 1) / 2; d++)
							{
								//std::cout << v[0] << "," << v[1] << "," << v[2] << "," << d << "\n";
								con1 = (d <= v[0]);
								con2 = (d <= v[1]);
								con3 = (v[0] + l < ns - d);
								con4 = (v[1] < ns - d);
								if (con1 && con2 && con3 && con4)
								{
									st_loops[l - 1][w - 1][d][config % arrsize] += wloop;
								}
							}
						}
					}
					for (v[0] = 0; v[0] < ns - 1; v[0]++)
					{
						for (v[1] = 0; v[1] < ns - l; v[1]++)
						{
							wloop = ((float)CalculateWilsonLoop(v, 1, 2, l, w));
							for (int d = 0; d < (ns - 1) / 2; d++)
							{
								con1 = (d <= v[0]);
								con2 = (d <= v[1]);
								con3 = (v[0] < ns - d );
								con4 = (v[1] + l < ns - d);
								if (con1 && con2 && con3 && con4)
								{
									st_loops[l - 1][w - 1][d][config % arrsize] += wloop;
								}
							}
						}
					}
				}
			}
		}
		if ((config + 1) % arrsize == 0)
		{

			std::ofstream outfile;
			std::string str1 = std::to_string(betatau);
			std::string str2 = std::to_string(coupling);
			str1 = str1.replace(str1.find('.'), 1, "_");
			str2 = str2.replace(str2.find('.'), 1, "_");
			for (int l = 0; l < (ns - 1) / 2; l++)
			{
				for (int w = 0; w < (ns - 1) / 2; w++)
				{
					filename = "../wilsonloops/";
					filename += "wilsonloops_ss_l=" + std::to_string(l) + "w=" + std::to_string(w) + "beta=" + str1;
					filename += "anisotropy=" + str2;
					filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
					filename += ".csv";
					std::cout << filename << "\n";
					outfile.open(filename, std::ios::app);
					for (int c = 0; c < arrsize; c++)
					{
						int d = 0;
						for (d = 0; d < (ns - 1) / 2 - 1; d++)
						{
							//outfile << ss_loops[l][w][d][c] / ((float)std::max(1, (ns - l - 1 - 2 * d) * (ns - w - 1 - 2 * d) * nt)) << ",";
							outfile << ss_loops[l][w][d][c] / ((float)ss_loops_norm[l][w][d]) << ",";
							ss_loops[l][w][d][c] = 0;
						}
						//outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)(ns - l - 2 * (d)) * (ns - 2 * d) * 2 * nt) << "\n";
						//outfile << ss_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)std::max(1, (ns - l - 1 - 2 * d)) * (ns - w - 1 - 2 * d) * nt) << "\n";
						outfile << ss_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)ss_loops_norm[l][w][(ns - 1) / 2 - 1]) << "\n";
						ss_loops[l][w][(ns - 1) / 2 - 1][c] = 0;
					}
					outfile.close();
				}
				std::cout << "spacial done\n";
				for (int w = 0; w < nt / 2; w++)
				{
					filename = "../wilsonloops/";
					filename += "wilsonloops_st_l=" + std::to_string(l) + "w=" + std::to_string(w) + "beta=" + str1;
					filename += "anisotropy=" + str2;
					filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
					filename += ".csv";
					outfile.open(filename, std::ios::app);
					for (int c = 0; c < arrsize; c++)
					{
						int d = 0;
						for (d = 0; d < (ns - 1) / 2 - 1; d++)
						{
							float norm = ((float)(ns - l - 1 - 2 * d) * (ns - 2 * d) * 2 * nt);
							//outfile << st_loops[l][w][d][c] / std::max(1.0f, norm) << ",";
							outfile << st_loops[l][w][d][c] / st_loops_norm[l][w][d] << ",";
							st_loops[l][w][d][c] = 0;
						}
						outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / st_loops_norm[l][w][(ns - 1) / 2 - 1] << "\n";
						//outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)std::max(1, (ns - l - 2 * d) * (ns - 2 * d) * 2 * nt)) << "\n";
						st_loops[l][w][(ns - 1) / 2 - 1][c] = 0;
					}
					outfile.close();
				}
			}
		}
	}
	std::ofstream outfile;
	std::string str1 = std::to_string(betatau);
	std::string str2 = std::to_string(coupling);
	str1 = str1.replace(str1.find('.'), 1, "_");
	str2 = str2.replace(str2.find('.'), 1, "_");
	for (int l = 0; l < (ns - 1) / 2; l++)
	{
		for (int w = 0; w < (ns - 1) / 2; w++)
		{
			filename = "../wilsonloops/";
			filename += "wilsonloops_ss_l=" + std::to_string(l) + "w=" + std::to_string(w) + "beta=" + str1;
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			filename += ".csv";
			outfile.open(filename, std::ios::app);
			for (int c = 0; c < configurations % arrsize; c++)
			{
				int d = 0;
				for (d = 0; d < (ns - 1) / 2 - 1; d++)
				{
					outfile << ss_loops[l][w][d][c] / ((float)(ns - l - 1 - 2 * d) * (ns - w - 1 - 2 * d) * nt)<< ",";
				}
				outfile << ss_loops[l][w][d][c] / ((float)(ns - l - 1 - 2 * d) * (ns - w - 1 - 2 * d) * nt) << "\n";
				//outfile << ss_loops[l][w][(ns - 1) / 2 - 1][c] / std::max(1.0f, ((float)(ns - l - 2 * (d + 1)) * (ns - w - 2 * (d + 1)) * nt)) << "\n";
			}
			outfile.close();
		}
		for (int w = 0; w < nt / 2; w++)
		{
			filename = "../wilsonloops/";
			filename += "wilsonloops_st_l=" + std::to_string(l) + "w=" + std::to_string(w) + "beta=" + str1;
			filename += "anisotropy=" + str2;
			filename += "ns_" + std::to_string(ns) + "_nt_" + std::to_string(nt);
			filename += ".csv";
			outfile.open(filename, std::ios::app);
			for (int c = 0; c < configurations % arrsize; c++)
			{
				int d = 0;
				for (d = 0; d < (ns - 1) / 2 - 1; d++)
				{
					float norm = ((float)(ns - l - 1 - 2 * d) * (ns - 2 * d) * 2 * nt);
					//outfile << st_loops[l][w][d][c] / std::max(1.0f, norm) << ",";
					outfile << st_loops[l][w][d][c] / st_loops_norm[l][w][d] << ",";
					st_loops[l][w][d][c] = 0;
				}
				outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / st_loops_norm[l][w][(ns - 1) / 2 - 1] << "\n";
				//outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)std::max(1, (ns - l - 2 * d) * (ns - 2 * d) * 2 * nt)) << "\n";
				st_loops[l][w][(ns - 1) / 2 - 1][c] = 0;
				//outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)(ns - l - 2 * (d)) * (ns - 2 * d) * 2 * nt) << "\n";
				//outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)(ns - l - 2 * d) * (ns - 2 * d) * 2 * nt) << "\n";
				//outfile << st_loops[l][w][(ns - 1) / 2 - 1][c] / ((float)std::max(1, (ns - l - 2 * d) * (ns - 2 * d) * 2 * nt)) << "\n";
			}
			outfile.close();
		}
	}
}

void Lattice::CalculateCorrelatorsOnTheFly(int ns, int nt, double coupling, double betatau, std::vector<std::vector<std::string>> operator_filenames, std::vector<std::vector<std::vector<std::string>>> correlatorfilenames)
{
	int x[3];
	int opernum = 4;
	std::cout << "Calculating operators\n";
	std::vector<std::vector<std::vector<double> > > operators;
	std::vector<std::vector<std::vector<double> > > correlators;
	operators.resize(opernum);
	correlators.resize(opernum);
	for (int i = 0; i < opernum; i++)
	{
		operators[i].resize(nt);
		for (int j = 0; j < nt; j++)
		{
			operators[i][j].resize((ns - 1) / 2);
		}
		correlators[i].resize(nt / 2);
		for (int t = 0; t < nt / 2; t++)
		{
			correlators[i][t].resize((ns - 1) / 2);
		}
	}
	for (int t = 0; t < nt; t++)
	{
		std::vector<std::vector<double> > operators_1 = CalculateOperators(t);
		for (int op1 = 0; op1 < opernum; op1++)
		{
			for (int d = 0; d < (ns - 1) / 2; d++)
			{
				operators[op1][t][d] += operators_1[op1][d];
			}
		}
	}
	for (int d = 0; d < (ns - 1) / 2; d++)
	{
		for (int t = 0; t < nt; t++)
		{
			for (int t0 = 0; t0 < nt / 2; t0++)
			{
				for (int op1 = 0; op1 < opernum; op1++)
				{
					correlators[op1][t0][d] += operators[op1][t][d] * operators[op1][(t + t0) % nt][d] / ((double)nt);
				}
			}
		}
	}
	std::ofstream outfile;
	for (int d = 0; d < (ns - 1) / 2; d++)
	{
		for (int op1 = 0; op1 < opernum; op1++)
		{
			float opval = 0.0;
			for (int t = 0; t < nt; t++)
			{
				opval += operators[op1][t][d] / ((float)nt);
			}
			outfile.open(operator_filenames[op1][d], std::ios::app);
			outfile << opval << "\n";
			outfile.close();\
			for (int dt = 0; dt < nt / 2; dt++)
			{
				outfile.open(correlatorfilenames[op1][dt][d], std::ios::app);
				outfile << correlators[op1][dt][d] << ", 0.0\n";
				outfile.close();
			}
		}
	}
}

int Lattice::CalculateWilsonLoop(int x[], int d, int dperp, int height, int width)
{
	/*x[] is the site,
	d is one direction of the wilson loop,
	dperp is a second direction of the wilson loop
	height is the length of the wilson loop in d direction
	width is the length of the  wilson loop in dperp direction*/
	// check if the wilson loop would go off the edge
	int wilsonloop = 1;
	for (int i = 0; i < height; i++)
	{
		wilsonloop *= link[GetIndex(x, d)];
		MoveUp(x, d);
	}
	for (int i = 0; i < width; i++) {
		wilsonloop *= link[GetIndex(x, dperp)];
		MoveUp(x, dperp);
	}
	for (int i = 0; i < height; i++) {
		MoveDown(x, d);
		wilsonloop *= link[GetIndex(x, d)];
	}
	for (int i = 0; i < width; i++) {
		MoveDown(x, dperp);
		wilsonloop *= link[GetIndex(x, dperp)];
	}
	//std::cout << wilsonloop << "\n";
	return wilsonloop;
}

void Lattice::CalculateAllSpacialLoops()
{
	int xshift = (ns - 1) / 4;
	int x[3];
	// iterate through the x sites away from 1/4 of the edges x side
	for (x[0] = xshift; x[0] < xshift + (ns - 1) / 2; x[0]++)
	{
		// iterate through the y sites away from 1/4 of the edges y side
		for (x[1] = xshift; x[1] < xshift + (ns - 1) / 2; x[1]++)
		{
			// iterate through the time steps
			for (x[2] = 0; x[2] < nt; x[2]++)
			{
				// iterate through the possible widths
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
							wilson_ss[w + (ns - 1) / 2 * l] += CalculateWilsonLoop(y, 0, 1, l + 1, w + 1) / norms_ss[w + (ns - 1) / 2 * l];
						}
						else
						{
							wilson_ss[l + (ns - 1) / 2 * w] += CalculateWilsonLoop(y, 0, 1, l + 1, w + 1) / norms_ss[l + (ns - 1) / 2 * w];
						}
					}
				}
			}
		}
	}
}

void Lattice::CalculateAllTemporalLoops()
{
	int xshift = (ns - 1) / 4;
	int x[3];
	// iterate through the x sites away from 1/4 of the edges x side
	for (x[0] = xshift; x[0] < xshift + (ns - 1) / 2; x[0]++)
	{
		// iterate through the y sites away from 1/4 of the edges y side
		for (x[1] = xshift; x[1] < xshift + (ns - 1) / 2; x[1]++)
		{
			// iterate through the time steps
			for (x[2] = 0; x[2] < nt; x[2]++)
			{
				// iterate through the possible widths
				for (int w = 0; w < (ns - 1) / 2 - (x[0] - xshift); w++)
				{
					for (int l = 0; l < nt / 2; l++)
					{
						int y[3];
						y[0] = 0 + x[0];
						y[1] = 0 + x[1];
						y[2] = 0 + x[2];
						wilson_st[w * (ns - 1) / 2 + l] += CalculateWilsonLoop(y, 0, 2, w, l) / norms_st[w * (ns - 1) / 2 + l];
					}
				}
				// iterate through the possible widths
				for (int w = 0; w < (ns - 1) / 2 - (x[1] - xshift); w++)
				{
					for (int l = 0; l < nt / 2; l++)
					{
						int y[3];
						y[0] = 0 + x[0];
						y[1] = 0 + x[1];
						y[2] = 0 + x[2];
						wilson_st[w * (ns - 1) / 2 + l] += CalculateWilsonLoop(y, 1, 2, w, l) / norms_st[w * (ns - 1) / 2 + l];
					}
				}
			}
		}
	}
}

void Lattice::CalculateWilsonLoopsOnTheFly(int ns, int nt, double coupling, double betatau, int xshift, std::vector<std::vector< std::string> > wilson_ss_filenames, std::vector<std::vector< std::string> > wilson_st_filenames, std::vector<std::string> correlator_filenames, bool calc_corr, bool calc_wil)
{
	if (calc_wil)
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
								loops_ss[w][l].push_back(CalculateWilsonLoop(y, 0, 1, l + 1, w + 1));
							}
							else
							{
								loops_ss[l][w].push_back(CalculateWilsonLoop(y, 0, 1, l + 1, w + 1));
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
							loops_st[w][l].push_back(CalculateWilsonLoop(y, 0, 2, l + 1, w + 1));
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
							loops_st[w][l].push_back(CalculateWilsonLoop(y, 1, 2, l + 1, w + 1));
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
	if (calc_corr)
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
						int y[3];
						int z[3];
						y[0] = 0 + x[0];
						y[1] = 0 + x[1];
						y[2] = 0 + x[2];
						z[0] = 0 + x[0];
						z[1] = 0 + x[1];
						z[2] = (w + x[2]) % nt;
						int t1 = CalculateWilsonLoop(x, 0, 1, 1, 1);
						int t2 = CalculateWilsonLoop(z, 0, 1, 1, 1);
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
	}
}

void Lattice::CalculatePlaquetteFromCenter(int ns, int nt, double coupling, double betatau, std::vector<std::vector<std::string> > plaqss_loc_filenames, std::vector<std::vector<std::vector<std::string> > >plaqst_loc_filenames)
{
	int x[3];
	int index;
	int plaq = 1;
	std::vector<std::vector<std::vector<int> > > plaq_ss;
	std::vector<std::vector<double> > mean_ss;
	std::vector<std::vector<std::vector<std::vector<int> > > > plaq_st;
	std::vector<std::vector<std::vector<double> > > mean_st;
	plaq_ss.resize(ns - 1);
	plaq_st.resize(ns);
	mean_ss.resize(ns - 1);
	mean_st.resize(ns);
	for (int i = 0; i < ns - 1; i++)
	{
		plaq_ss[i].resize(ns - 1);
		mean_ss[i].resize(ns - 1);
	}
	for (int i = 0; i < ns; i++)
	{
		plaq_st[i].resize(ns);
		mean_st[i].resize(ns);
		for (int j = 0; j < ns; j++)
		{
			plaq_st[i][j].resize(2);
			mean_st[i][j].resize(2);
		}
	}
	for (x[0] = 0; x[0] < ns; x[0]++)
	{
		for (x[1] = 0; x[1] < ns; x[1]++)
		{
			for (x[2] = 0; x[2] < ns; x[2]++)
			{
				plaq = 1;
				// only if not at the right end of lattice
				if (x[1] < ns - 1 && x[0] < ns - 1)
				{
					index = GetIndex(x, 0);
					plaq *= link[index];
					MoveUp(x, 0);
					plaq *= link[GetIndex(x, 1)];
					MoveUp(x, 1);
					MoveDown(x, 0);
					plaq *= link[GetIndex(x, 0)];
					MoveDown(x, 1);
					plaq *= link[GetIndex(x, 1)];
					plaq_ss[x[0]][x[1]].push_back(0 + plaq);
				}
				plaq = 1;
				if (x[0] < ns - 1)
				{
					plaq *= link[GetIndex(x, 0)];
					MoveUp(x, 0);
					plaq *= link[GetIndex(x, 2)];
					MoveUp(x, 2);
					MoveDown(x, 0);
					plaq *= link[GetIndex(x, 0)];
					MoveDown(x, 2);
					plaq *= link[GetIndex(x, 2)];
					plaq_st[x[0]][x[1]][0].push_back(plaq + 0);
				}
				plaq = 1;
				if (x[1] < ns - 1)
				{
					plaq *= link[GetIndex(x, 1)];
					MoveUp(x, 1);
					plaq *= link[GetIndex(x, 2)];
					MoveUp(x, 2);
					MoveDown(x, 1);
					plaq *= link[GetIndex(x, 1)];
					MoveDown(x, 2);
					plaq *= link[GetIndex(x, 2)];
					plaq_st[x[0]][x[1]][1].push_back(plaq + 0);
				}
			}
		}
	}
	for (x[0] = 0; x[0] < ns - 1; x[0]++)
	{
		for (x[1] = 0; x[1] < ns - 1; x[1]++)
		{
			double num_elements = ((double)plaq_ss[x[0]][x[1]].size());
			for (int i = 0; i < num_elements; i++)
			{
				mean_ss[x[0]][x[1]] += ((double)plaq_ss[x[0]][x[1]][i]) / num_elements;
			}
		}
	}
	for (x[0] = 0; x[0] < ns; x[0]++)
	{
		for (x[1] = 0; x[1] < ns; x[1]++)
		{
			if (x[0] < ns - 1)
			{
				double num_elements = ((double)plaq_st[x[0]][x[1]][0].size());
				for (int i = 0; i < num_elements; i++)
				{
					mean_st[x[0]][x[1]][0] += ((double)plaq_st[x[0]][x[1]][0][i]) / num_elements;
				}
			}
			if (x[1] < ns - 1)
			{
				double num_elements = ((double)plaq_st[x[0]][x[1]][1].size());
				for (int i = 0; i < num_elements; i++)
				{
					mean_st[x[0]][x[1]][1] += ((double)plaq_st[x[0]][x[1]][1][i]) / num_elements;
				}
			}
		}
	}
	std::ofstream outfile;
	for (x[0] = 0; x[0] < ns - 1; x[0]++)
	{
		for (x[1] = 0; x[1] < ns - 1; x[1]++)
		{
			outfile.open(plaqss_loc_filenames[x[0]][x[1]], std::ios::app);
			outfile << mean_ss[x[0]][x[1]] << "," << 0 << "\n";
			outfile.close();
		}
	}
	for (x[0] = 0; x[0] < ns; x[0]++)
	{
		for (x[1] = 0; x[1] < ns; x[1]++)
		{
			if (x[0] < ns - 1)
			{
				outfile.open(plaqst_loc_filenames[x[0]][x[1]][0], std::ios::app);
				outfile << mean_st[x[0]][x[1]][0] << "," << 0 << "\n";
				outfile.close();
			}
			if (x[1] < ns - 1)
			{
				outfile.open(plaqst_loc_filenames[x[0]][x[1]][1], std::ios::app);
				outfile << mean_st[x[0]][x[1]][1] << "," << 0 << "\n";
				outfile.close();
			}
		}
	}
}
