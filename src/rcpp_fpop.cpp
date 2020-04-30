#include <Rcpp.h>
#include "Fpop.h"
#include <string>
#include <vector>
#include <map>
using namespace Rcpp;


// [[Rcpp::export]]
List fpop_cpp(std::vector<double> y, double alpha, double muMinLocal, double muMaxLocal, std::vector<double>  wt)
{   
  Fpop f = Fpop (y, 
    alpha, 
    muMinLocal,
    muMaxLocal,
    wt
  );
      
  f.Search();
      
  List l = List::create(
    _["changepoints"] = f.Retreive_changepoints(),
    _["costs"] = f.Retreive_costs(),
    _["means"] = f.Retreive_means(),
    _["intervals"] = f.Get_intervals(),
    _["candidates"] = f.Get_candidates()
  );
      
  return l;
}
