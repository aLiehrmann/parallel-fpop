#include <Rcpp.h>
#include "Fpop.h"
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include "omp.h"

using namespace Rcpp;

// [[Rcpp::export]]
#ifndef CPP
List fpop_cpp(std::vector<double> y, double alpha, double muMinLocal, double muMaxLocal , std::vector<double> wt, int nbThreads, int algorithmChoice)
#else
void fpop_cpp(std::vector<double> y, double alpha, double muMinLocal, double muMaxLocal , std::vector<double> wt, int nbThreads, int algorithmChoice) //, std::vector<double> v)
#endif
{
    // List l;

    if (1 == nbThreads)
    {
        // Sequential version of the algorithm, do not consider the algorithmChoice parameter
        double tStart = omp_get_wtime();
        Fpop f = Fpop(y, alpha, muMinLocal, muMaxLocal, wt);
        f.Search();
        double tEnd = omp_get_wtime();

        #ifndef CPP
            l.push_back(List::create(_["id"] = 0,
                _["changepoints"] = f.Retrieve_changepoints(),
                _["costs"] = f.Retrieve_costs(),
                _["means"] = f.Retrieve_means(),
                _["intervals"] = f.Get_intervals(),
                _["candidates"] = f.Get_candidates(),
                _["execution_time"] = tEnd - tStart
            ));
        #else
            printf("Time: %f", tEnd - tStart);

            std::vector<int> changepoints = f.Retrieve_changepoints();
            std::vector<double> costs = f.Retrieve_costs();
            std::vector<double> means = f.Retrieve_means();
            std::vector<int> candidates = f.Get_candidates();
            std::vector<int> intervals = f.Get_intervals();

            printf("\nChangepoints: ");
            for (int i = 0; i < changepoints.size(); ++i)
            {
                printf("%d ", changepoints[i]);
            }
            printf("\nCosts: ");
            for (int i = 0; i < costs.size(); ++i)
            {
                printf("%f ", costs[i]);
            }
            printf("\nMeans: ");
            for (int i = 0; i < means.size(); ++i)
            {
                printf("%f ", means[i]);
            }
            printf("\nCandidates: %d", candidates.size());
            printf("\nIntervals: %d", intervals.size());
        #endif

    } else {
        // Parallel version of the algorithm
        if (0 == algorithmChoice || 1 == algorithmChoice)
        {
            // Declare variables that are common to whether we partition the work based on strict
            // partition of mu or based on the number of candidate points
            Fpop * f = new Fpop[nbThreads];
            double * F = new double[nbThreads];
            double * ARG_F = new double[nbThreads];
            int * t_hat = new int[nbThreads];

            // Variables related to "reducing" the local minimum point of each thread to the global minimum point
            bool * firstMin = new bool;
            (*firstMin) = false;
            double * F_min = new double;
            (*F_min) = 0.0;
            double * ARG_F_min = new double;
            (*ARG_F_min) = 0.0;
            int * t_hat_min = new int;
            (*t_hat_min) = 0;

            if (0 == algorithmChoice)
            {
                // The work is uniformly split across mu, i.e., mu is split in as many intervals as there are threads
                double muRange = (muMaxLocal - muMinLocal) / nbThreads;

                #pragma omp parallel num_threads(nbThreads)
                {
                    double tStart = omp_get_wtime();
                    int tid = omp_get_thread_num();

                    double muMinLocal_ = muMinLocal + muRange * tid;
                    double muMaxLocal_ = muMinLocal_ + muRange;

                    f[tid] = Fpop(y, alpha, muMinLocal_, muMaxLocal_, wt);
                    //f[tid] = Fpop(y, alpha, v[tid], v[tid+1], wt);
                    f[tid].Search_parallel(tid, nbThreads, F, ARG_F, t_hat, firstMin, F_min, ARG_F_min, t_hat_min);
                    double tEnd = omp_get_wtime();

                    #ifndef CPP
                        #pragma omp critical
                        l.push_back(List::create(_["id"] = tid,
                            _["changepoints"] = f[tid].Retrieve_changepoints(),
                            _["costs"] = f[tid].Retrieve_costs(),
                            _["means"] = f[tid].Retrieve_means(),
                            _["intervals"] = f[tid].Get_intervals(),
                            _["candidates"] = f[tid].Get_candidates(),
                            _["execution_time"] = tEnd - tStart
                        ));
                    #else
                        printf("Time: %f\n", tEnd - tStart);
                    #endif
                }
            } else {
                // The work is split across the threads based on the number of candidate points.
                // I.e., partition mu to attribute a uniform number of candidate points to each thread
                int nbPoints = y.size();
                int pointsPerThread = nbPoints / nbThreads;
                int threadsWithOneMore = nbPoints - (pointsPerThread * nbThreads);

                #pragma omp parallel num_threads(nbThreads)
                {
                    double tStart = omp_get_wtime();
                    int tid = omp_get_thread_num();

                    // Offset the real beginning and ending of the interval considering the threads
                    // with one more point (to comply with uneven divisions)
                    int offsetMin;
                    int offsetMax;
                    // Ensure that we are not going out of bounds
                    if (tid < threadsWithOneMore)
                    {
                        // Tid 0 offsets 0, tid 1 offset 1 (from 0), etc.
                        offsetMin = tid;
                        offsetMax = 1;
                    } else {
                        // Offset all the threads with one more point in the interval
                        offsetMin = threadsWithOneMore;
                        offsetMax = 0;
                    }
                    // Get the boundaries of the interval
                    double muMin = y[tid * pointsPerThread + offsetMin];
                    double muMax = y[tid * pointsPerThread + pointsPerThread + offsetMax];

                    printf("Interval %d: [%f, %f]\n", tid, muMin, muMax);

                    f[tid] = Fpop(y, alpha, muMin, muMax, wt);
                    //f[tid] = Fpop(y, alpha, v[tid], v[tid+1], wt);
                    f[tid].Search_parallel(tid, nbThreads, F, ARG_F, t_hat, firstMin, F_min, ARG_F_min, t_hat_min);
                    double tEnd = omp_get_wtime();

                    #ifndef CPP
                        #pragma omp critical
                        l.push_back(List::create(_["id"] = tid,
                            _["changepoints"] = f[tid].Retrieve_changepoints(),
                            _["costs"] = f[tid].Retrieve_costs(),
                            _["means"] = f[tid].Retrieve_means(),
                            _["intervals"] = f[tid].Get_intervals(),
                            _["candidates"] = f[tid].Get_candidates(),
                            _["execution_time"] = tEnd - tStart
                        ));
                    #else
                        printf("Time: %f\n", tEnd - tStart);
                    #endif
                }
            }

            #ifdef CPP
                for (int i = 0; i < nbThreads; ++i)
                {
                    std::vector<int> changepoints = f[i].Retrieve_changepoints();
                    std::vector<double> costs = f[i].Retrieve_costs();
                    std::vector<double> means = f[i].Retrieve_means();
                    std::vector<int> candidates = f[i].Get_candidates();
                    std::vector<int> intervals = f[i].Get_intervals();

                    printf("Thread %d: \n", i);

                    printf("Changepoints: ");
                    for (int i = 0; i < changepoints.size(); ++i)
                    {
                        printf("%d ", changepoints[i]);
                    }
                    printf("\nCosts: ");
                    for (int i = 0; i < costs.size(); ++i)
                    {
                        printf("%f ", costs[i]);
                    }
                    printf("\nMeans: ");
                    for (int i = 0; i < means.size(); ++i)
                    {
                        printf("%f ", means[i]);
                    }
                    printf("\nCandidates: %d", candidates.size());
                    printf("\nIntervals: %d", intervals.size());
                    printf("\n\n");
                }
            #endif

            delete[] f;
            delete[] F;
            delete[] ARG_F;
            delete[] t_hat;

            delete firstMin;
            delete F_min;
            delete ARG_F_min;
            delete t_hat_min;

        } else {
            if (2 == algorithmChoice)
            {
                // Sequential main loop to find the minimizing point, with several parallelized loops
                double tStart = omp_get_wtime();
                Fpop f = Fpop(y, alpha, muMinLocal, muMaxLocal, wt);
                f.Search_parallel_loops(nbThreads);
                double tEnd = omp_get_wtime();

                #ifndef CPP
                    l.push_back(List::create(_["id"] = 0,
                        _["changepoints"] = f.Retrieve_changepoints(),
                        _["costs"] = f.Retrieve_costs(),
                        _["means"] = f.Retrieve_means(),
                        _["intervals"] = f.Get_intervals(),
                        _["candidates"] = f.Get_candidates(),
                        _["execution_time"] = tEnd - tStart
                    ));
                #else
                    printf("Time: %f", tEnd - tStart);

                    std::vector<int> changepoints = f.Retrieve_changepoints();
                    std::vector<double> costs = f.Retrieve_costs();
                    std::vector<double> means = f.Retrieve_means();
                    std::vector<int> candidates = f.Get_candidates();
                    std::vector<int> intervals = f.Get_intervals();

                    printf("\nChangepoints: ");
                    for (int i = 0; i < changepoints.size(); ++i)
                    {
                        printf("%d ", changepoints[i]);
                    }
                    printf("\nCosts: ");
                    for (int i = 0; i < costs.size(); ++i)
                    {
                        printf("%f ", costs[i]);
                    }
                    printf("\nMeans: ");
                    for (int i = 0; i < means.size(); ++i)
                    {
                        printf("%f ", means[i]);
                    }
                    printf("\nCandidates: %d", candidates.size());
                    printf("\nIntervals: %d", intervals.size());
                #endif
            } else {
                std::cout << "Error : Unknown algorithm choice\n";
            }
        }
    }

    // return l;
}
