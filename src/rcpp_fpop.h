#ifndef RCPP_FPOP_H
#define RCPP_FPOP_H

#ifndef R
void fpop_cpp(std::vector<double> y, double alpha, double muMinLocal, double muMaxLocal , std::vector<double> wt, int nbThreads, int algorithmChoice);
#else
List fpop_cpp(std::vector<double> y, double alpha, double muMinLocal, double muMaxLocal , std::vector<double> wt, int nbThreads, int algorithmChoice);
#endif

#endif
