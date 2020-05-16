#include "Fpop.h"
#include "rcpp_fpop.h"

#include <random>
#include <iostream>
#include <vector>

#include "omp.h"

int main(int argc, char * argv[])
{
    if (argc != 4)
    {
        std::cout << "Incorrect number of parameters (expected: 3)\n";
        return 1;
    }

    int size = atoi(argv[1]);
    int nbThreads = atoi(argv[2]);
    int algorithmChoice = atoi(argv[3]);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);

    std::vector<double> vectData;
    for (int i = 0; i < size; ++i)
    {
        vectData.push_back(distribution(generator));
    }

    std::vector<double> weight(size, 1.0);

    double muMin = 0.0;
    double muMax = 0.0;

    for (int i = 0; i < size; ++i)
    {
        if (vectData[i] < muMin)
        {
            muMin = vectData[i];
        }
        if (muMax < vectData[i])
        {
            muMax = vectData[i];
        }
    }

    printf("Mu min = %f, mu max = %f\n", muMin, muMax);

    fpop_cpp(vectData, 2.5, muMin, muMax, weight, nbThreads, algorithmChoice);
}
