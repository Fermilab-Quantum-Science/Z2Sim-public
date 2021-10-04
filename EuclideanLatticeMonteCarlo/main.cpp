#include <iostream>
#include <stdlib.h>
#include "Lattice.h"
#include "configgeneration.h"
#include "WilsonLoopCalc.h"

/******************************/
int main(int argc, char* argvp[])
{
    std::cout << argvp[1] << "\n";
    if (strcmp(argvp[1], "-h") == 0)
    {
        std::cout << "print help line late";
    }
    // if first flag is 'gen'  then generate configurations
    if (strcmp(argvp[1], "gen") == 0)
    {
        if (argc < 5)
        {
            std::cout << "insufficient number of arguments\n";
        }
        else
        {
            int ns = atoi(argvp[2]);
            int nt = atoi(argvp[3]);
            double anisotropy = atof(argvp[4]);
            double beta = atof(argvp[5]);
            int burninsteps = 1000;
            int montecarlosweeps = 100;
            int configurations = 1;
            int startconfig = 0;
            bool calcwilson = false;
            for (int i = 5; i < argc; i++)
            {
                if (strcmp(argvp[i], "bis") == 0)
                {
                    burninsteps = atoi(argvp[i + 1]);
                }
                if (strcmp(argvp[i], "mcs") == 0)
                {
                    montecarlosweeps = atoi(argvp[i + 1]);
                }
                if (strcmp(argvp[i], "ncon") == 0)
                {
                    configurations = atoi(argvp[i + 1]);
                }
                if (strcmp(argvp[i], "start") == 0)
                {
                    startconfig = atoi(argvp[i + 1]);
                }
                if (strcmp(argvp[i], "wilson") == 0)
                {
                    calcwilson = false;
                }
            }
            std::cout << "generating configurations:\n";
            std::cout << "\t ns: \t" << ns << "\n";
            std::cout << "\t nt: \t" << nt << "\n";
            std::cout << "\t bare anisotropy: \t" << anisotropy << "\n";
            std::cout << "\t beta: \t" << beta << "\n";
            std::cout << "\t burn in steps: \t" << burninsteps << "\n";
            std::cout << "\t monte carlo sweeps: \t" << montecarlosweeps << "\n";
            std::cout << "\t configurations to generate: \t" << configurations << "\n";
            srand(time(NULL));
            generate_configurations(ns, nt, burninsteps, montecarlosweeps, configurations, anisotropy, beta, startconfig);
        }
    }
    // check if we want to generate wilson loops
    if (strcmp(argvp[1], "wilson") == 0)
    {
        if (argc < 5)
        {
            std::cout << "insufficient number of arguments\n";
        }
        else
        {
            int ns = atoi(argvp[2]);
            int nt = atoi(argvp[3]);
            double magcoup = atof(argvp[4]);
            double betatau = atof(argvp[5]);
            int configurationsstart = 0;
            int configurationsend = 100;
            int set = 0;
            for (int i = 5; i < argc; i++)
            {
                if (strcmp(argvp[i], "start") == 0)
                {
                    configurationsstart = atoi(argvp[i + 1]);
                }
                if (strcmp(argvp[i], "stop") == 0)
                {
                    configurationsend = atoi(argvp[i + 1]);
                }
                if (strcmp(argvp[i], "set") == 0)
                {
                    set = atoi(argvp[i + 1]);
                }
            }
            std::cout << "calculating Wilson loops:\n";
            std::cout << "\t ns: \t" << ns << "\n";
            std::cout << "\t nt: \t" << nt << "\n";
            std::cout << "\t anisotropy: \t" << magcoup << "\n";
            std::cout << "\t beta: \t" << betatau << "\n";
            std::cout << "\t for configurations: \t" << configurationsstart << " to " << configurationsend << "\n";
            CalculateWilsonLoopsFromConfig(ns, nt, configurationsstart, configurationsend, magcoup, betatau, set);
        }
    }
    if (strcmp(argvp[1], "corr") == 0)
    {
        if (argc < 5)
        {
            std::cout << "insufficient number of arguments\n";
        }
        else
        {
            int ns = atoi(argvp[2]);
            int nt = atoi(argvp[3]);
            double magcoup = atof(argvp[4]);
            double betatau = atof(argvp[5]);
            int configurationsstart = 0;
            int configurationsend = 100;
            int set = 0;
            for (int i = 5; i < argc; i++)
            {
                if (strcmp(argvp[i], "start") == 0)
                {
                    configurationsstart = atoi(argvp[i + 1]);
                }
                if (strcmp(argvp[i], "stop") == 0)
                {
                    configurationsend = atoi(argvp[i + 1]);
                }
                if (strcmp(argvp[i], "set") == 0)
                {
                    set = atoi(argvp[i + 1]);
                }
            }
            std::cout << "calculating correlators:\n";
            std::cout << "\t ns: \t" << ns << "\n";
            std::cout << "\t nt: \t" << nt << "\n";
            std::cout << "\t anisotropy: \t" << magcoup << "\n";
            std::cout << "\t beta: \t" << betatau << "\n";
            std::cout << "\t for configurations: \t" << configurationsstart << " to " << configurationsend << "\n";
            MeasureCorrelatorFromConfig(ns, nt, configurationsstart, configurationsend, magcoup, betatau, set);
        }
    }
    if (strcmp(argvp[1], "gw") == 0)
    {
        int ns = atoi(argvp[2]);
        int nt = atoi(argvp[3]);
        double anisotropy = atof(argvp[4]);
        double beta = atof(argvp[5]);
        int burninsteps = 1000;
        int montecarlosweeps = 100;
        int configurations = 1;
        int startconfig = 0;
        bool calcwilson = false;
        for (int i = 5; i < argc; i++)
        {
            if (strcmp(argvp[i], "bis") == 0)
            {
                burninsteps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "mcs") == 0)
            {
                montecarlosweeps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "ncon") == 0)
            {
                configurations = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "start") == 0)
            {
                startconfig = atoi(argvp[i + 1]);
            }
        }
        std::cout << "generating configurations:\n";
        std::cout << "\t ns: \t" << ns << "\n";
        std::cout << "\t nt: \t" << nt << "\n";
        std::cout << "\t bare anisotropy: \t" << anisotropy << "\n";
        std::cout << "\t beta: \t" << beta << "\n";
        std::cout << "\t burn in steps: \t" << burninsteps << "\n";
        std::cout << "\t monte carlo sweeps: \t" << montecarlosweeps << "\n";
        std::cout << "\t configurations to generate: \t" << configurations << "\n";
        srand(time(NULL));
        generate_configurations_and_calc(ns, nt, burninsteps, montecarlosweeps, configurations, anisotropy, beta, (ns - 1) / 4, true, false, true);
    }
    if (strcmp(argvp[1], "gcw") == 0)
    {
        int ns = atoi(argvp[2]);
        int nt = atoi(argvp[3]);
        double anisotropy = atof(argvp[4]);
        double beta = atof(argvp[5]);
        int burninsteps = 1000;
        int montecarlosweeps = 100;
        int configurations = 1;
        int startconfig = 0;
        bool calcwilson = false;
        for (int i = 5; i < argc; i++)
        {
            if (strcmp(argvp[i], "bis") == 0)
            {
                burninsteps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "mcs") == 0)
            {
                montecarlosweeps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "ncon") == 0)
            {
                configurations = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "start") == 0)
            {
                startconfig = atoi(argvp[i + 1]);
            }
        }
        std::cout << "generating configurations:\n";
        std::cout << "\t ns: \t" << ns << "\n";
        std::cout << "\t nt: \t" << nt << "\n";
        std::cout << "\t bare anisotropy: \t" << anisotropy << "\n";
        std::cout << "\t beta: \t" << beta << "\n";
        std::cout << "\t burn in steps: \t" << burninsteps << "\n";
        std::cout << "\t monte carlo sweeps: \t" << montecarlosweeps << "\n";
        std::cout << "\t configurations to generate: \t" << configurations << "\n";
        srand(time(NULL));
        generate_configurations_and_calc(ns, nt, burninsteps, montecarlosweeps, configurations, anisotropy, beta, (ns - 1) / 4, true, true, true);
    }
    if (strcmp(argvp[1], "gwl") == 0)
    {
        int ns = atoi(argvp[2]);
        int nt = atoi(argvp[3]);
        double anisotropy = atof(argvp[4]);
        double beta = atof(argvp[5]);
        int burninsteps = 1000;
        int montecarlosweeps = 100;
        int configurations = 1;
        int startconfig = 0;
        bool calcwilson = false;
        for (int i = 5; i < argc; i++)
        {
            if (strcmp(argvp[i], "bis") == 0)
            {
                burninsteps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "mcs") == 0)
            {
                montecarlosweeps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "ncon") == 0)
            {
                configurations = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "start") == 0)
            {
                startconfig = atoi(argvp[i + 1]);
            }
        }
        std::cout << "generating configurations:\n";
        std::cout << "\t ns: \t" << ns << "\n";
        std::cout << "\t nt: \t" << nt << "\n";
        std::cout << "\t bare anisotropy: \t" << anisotropy << "\n";
        std::cout << "\t beta: \t" << beta << "\n";
        std::cout << "\t burn in steps: \t" << burninsteps << "\n";
        std::cout << "\t monte carlo sweeps: \t" << montecarlosweeps << "\n";
        std::cout << "\t configurations to generate: \t" << configurations << "\n";
        srand(time(NULL));
        generate_configurations_and_calc_plaq(ns, nt, burninsteps, montecarlosweeps, configurations, anisotropy, beta, 0, true);
    }
    if (strcmp(argvp[1], "gcl") == 0)
    {
        int ns = atoi(argvp[2]);
        int nt = atoi(argvp[3]);
        double anisotropy = atof(argvp[4]);
        double beta = atof(argvp[5]);
        int burninsteps = 1000;
        int montecarlosweeps = 100;
        int configurations = 1;
        int startconfig = 0;
        bool calcwilson = false;
        for (int i = 5; i < argc; i++)
        {
            if (strcmp(argvp[i], "bis") == 0)
            {
                burninsteps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "mcs") == 0)
            {
                montecarlosweeps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "ncon") == 0)
            {
                configurations = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "start") == 0)
            {
                startconfig = atoi(argvp[i + 1]);
            }
        }
        std::cout << "generating configurations:\n";
        std::cout << "\t ns: \t" << ns << "\n";
        std::cout << "\t nt: \t" << nt << "\n";
        std::cout << "\t bare anisotropy: \t" << anisotropy << "\n";
        std::cout << "\t beta: \t" << beta << "\n";
        std::cout << "\t burn in steps: \t" << burninsteps << "\n";
        std::cout << "\t monte carlo sweeps: \t" << montecarlosweeps << "\n";
        std::cout << "\t configurations to generate: \t" << configurations << "\n";
        srand(time(NULL));
        generate_configurations_and_calc_cross_corr(ns, nt, burninsteps, montecarlosweeps, configurations, anisotropy, beta);
    }
    if (strcmp(argvp[1], "gdc") == 0)
    {
        int ns = atoi(argvp[2]);
        int nt = atoi(argvp[3]);
        double anisotropy = atof(argvp[4]);
        double beta = atof(argvp[5]);
        int burninsteps = 1000;
        int montecarlosweeps = 100;
        int configurations = 1;
        int startconfig = 0;
        bool calcwilson = false;
        for (int i = 5; i < argc; i++)
        {
            if (strcmp(argvp[i], "bis") == 0)
            {
                burninsteps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "mcs") == 0)
            {
                montecarlosweeps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "ncon") == 0)
            {
                configurations = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "start") == 0)
            {
                startconfig = atoi(argvp[i + 1]);
            }
        }
        std::cout << "generating configurations:\n";
        std::cout << "\t ns: \t" << ns << "\n";
        std::cout << "\t nt: \t" << nt << "\n";
        std::cout << "\t bare anisotropy: \t" << anisotropy << "\n";
        std::cout << "\t beta: \t" << beta << "\n";
        std::cout << "\t burn in steps: \t" << burninsteps << "\n";
        std::cout << "\t monte carlo sweeps: \t" << montecarlosweeps << "\n";
        std::cout << "\t configurations to generate: \t" << configurations << "\n";
        srand(time(NULL));
        generate_configurations_and_calc_corr(ns, nt, burninsteps, montecarlosweeps, configurations, anisotropy, beta);
    }
    if (strcmp(argvp[1], "wcloc") == 0)
    {
        int ns = atoi(argvp[2]);
        int nt = atoi(argvp[3]);
        double anisotropy = atof(argvp[4]);
        double beta = atof(argvp[5]);
        int burninsteps = 1000;
        int montecarlosweeps = 100;
        int configurations = 1;
        int startconfig = 0;
        bool calcwilson = false;
        for (int i = 5; i < argc; i++)
        {
            if (strcmp(argvp[i], "bis") == 0)
            {
                burninsteps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "mcs") == 0)
            {
                montecarlosweeps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "ncon") == 0)
            {
                configurations = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "start") == 0)
            {
                startconfig = atoi(argvp[i + 1]);
            }
        }
        std::cout << "generating configurations:\n";
        std::cout << "\t ns: \t" << ns << "\n";
        std::cout << "\t nt: \t" << nt << "\n";
        std::cout << "\t bare anisotropy: \t" << anisotropy << "\n";
        std::cout << "\t beta: \t" << beta << "\n";
        std::cout << "\t burn in steps: \t" << burninsteps << "\n";
        std::cout << "\t monte carlo sweeps: \t" << montecarlosweeps << "\n";
        std::cout << "\t configurations to generate: \t" << configurations << "\n";
        srand(time(NULL));
        generate_configurations_and_calc_spacial_wilson_and_corr(ns, nt, burninsteps, montecarlosweeps, configurations, anisotropy, beta);
    }
    if (strcmp(argvp[1], "wloc") == 0)
    {
        int ns = atoi(argvp[2]);
        int nt = atoi(argvp[3]);
        double anisotropy = atof(argvp[4]);
        double beta = atof(argvp[5]);
        int burninsteps = 1000;
        int montecarlosweeps = 100;
        int configurations = 1;
        int startconfig = 0;
        bool calcwilson = false;
        for (int i = 5; i < argc; i++)
        {
            if (strcmp(argvp[i], "bis") == 0)
            {
                burninsteps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "mcs") == 0)
            {
                montecarlosweeps = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "ncon") == 0)
            {
                configurations = atoi(argvp[i + 1]);
            }
            if (strcmp(argvp[i], "start") == 0)
            {
                startconfig = atoi(argvp[i + 1]);
            }
        }
        std::cout << "generating configurations:\n";
        std::cout << "\t ns: \t" << ns << "\n";
        std::cout << "\t nt: \t" << nt << "\n";
        std::cout << "\t bare anisotropy: \t" << anisotropy << "\n";
        std::cout << "\t beta: \t" << beta << "\n";
        std::cout << "\t burn in steps: \t" << burninsteps << "\n";
        std::cout << "\t monte carlo sweeps: \t" << montecarlosweeps << "\n";
        std::cout << "\t configurations to generate: \t" << configurations << "\n";
        srand(time(NULL));
        generate_configurations_and_calc_spacial_wilson(ns, nt, burninsteps, montecarlosweeps, configurations, anisotropy, beta);
    }
}