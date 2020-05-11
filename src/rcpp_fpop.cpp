#include <Rcpp.h>
#include "Fpop.h"
#include <string>
#include <vector>
#include <map>
#include "omp.h"

using namespace Rcpp;

// [[Rcpp::export]]
List fpop_cpp(std::vector<double> y, double alpha, double muMinLocal, double muMaxLocal , std::vector<double> wt, int nbThreads) //, std::vector<double> v)
{
    Fpop f[nbThreads];
    double muRange = (muMaxLocal - muMinLocal)/nbThreads;
    double * F = new double[nbThreads];
    double * ARG_F = new double[nbThreads];
    int * t_hat = new int[nbThreads];

    double firstMin = false;
    double F_min = 0.0;
    double ARG_F_min = 0.0;
    int t_hat_min = 0;

    List l;
    #pragma omp parallel num_threads(nbThreads)
    {
        double tStart = omp_get_wtime();
        int tid = omp_get_thread_num();
        double muMinLocal_ = muMinLocal + muRange * tid;
        double muMaxLocal_ = muMinLocal_ + muRange;
        f[tid] = Fpop(y, alpha, muMinLocal_, muMaxLocal_, wt);
        //f[tid] = Fpop(y, alpha, v[tid], v[tid+1], wt);
        f[tid].Search(tid, nbThreads, F, ARG_F, t_hat, &firstMin, &F_min, &ARG_F_min, &t_hat_min);
        double tEnd = omp_get_wtime();

        #pragma omp critical
        l.push_back(List::create(_["id"] = tid,
            _["changepoints"] = f[tid].Retrieve_changepoints(),
            _["costs"] = f[tid].Retrieve_costs(),
            _["means"] = f[tid].Retrieve_means(),
            _["intervals"] = f[tid].Get_intervals(),
            _["candidates"] = f[tid].Get_candidates(),
            _["execution_time"] = tEnd - tStart
        ));
    }

    delete[] F;
    delete[] ARG_F;
    delete[] t_hat;

    return l;
}
