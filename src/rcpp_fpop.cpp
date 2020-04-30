#include <Rcpp.h>
#include "Fpop.h"
#include <string>
#include <vector>
#include <map>

#include "omp.h"

using namespace Rcpp;

// [[Rcpp::export]]
List * fpop_cpp(std::vector<double> y, double alpha, double muMinLocal, double muMaxLocal, std::vector<double> wt, int nbThreads)
{
    Fpop f[nbThreads];
    // Fpop f = Fpop (y, alpha, muMinLocal, muMaxLocal, wt);

    double muRange = (std::abs(muMinLocal) + std::abs(muMinLocal))/2;

    double * F = new double[nbThreads];
    double * ARG_F = new double[nbThreads];
    int * t_hat = new int[nbThreads];

    List * l[nbThreads];

    #pragma omp parallel num_threads(nbThreads)
    {
        int tid = omp_get_thread_num();

        double muMinLocal = muMin + muRange * tid;
        double muMaxLocal = muMinLocal + muRange;

        f[tid] = Fpop(y, alpha, muMinLocal, muMaxLocal, wt);

        f[tid].Search(tid, nbThreads, F, ARG_F, t_hat);

        l[tid] = List::create(
            _["changepoints"] = f[tid].Retrieve_changepoints(),
            _["costs"] = f[tid].Retrieve_costs(),
            _["means"] = f[tid].Retrieve_means(),
            _["intervals"] = f[tid].Get_intervals(),
            _["candidates"] = f[tid].Get_candidates()
        )
    }

    // List l = List::create(
    //     _["changepoints"] = f.Retrieve_changepoints(),
    //     _["costs"] = f.Retrieve_costs(),
    //     _["means"] = f.Retrieve_means(),
    //     _["intervals"] = f.Get_intervals(),
    //     _["candidates"] = f.Get_candidates()
    // );

    delete[] F;
    delete[] ARG_F;
    delete[] t_hat;

    return l;
}
