#ifndef DEF_FPOP
#define DEF_FPOP

#include <vector>
#include "Interval.h"
#include <random>

class Fpop
{
private:
    std::vector<double> y;
    int n; 
    Interval d;
    double lambda;
    std::vector<double> wt;
    std::vector<int> cp;
    std::vector<double> costs;
    std::vector<double> means;
    std::vector<int> nb_candidates;
    std::vector<int> nb_intervals;

    
public:
    /**
     * @param[in] y a data vector ordered according to an attribute.
     * @param[in] lambda a constant used in the calculation of the penalty.
     * @param[in] wt a vector of weight linked to the data.    
     */
    Fpop(std::vector<double> y_,
        double lambda_,
        double muMinLocal,
        double muMaxLocal, 
        std::vector<double> wt_ = {0});
    

    /**
     * @details Procedure for inferring the number of changepoints and their location in the data.
     */
    void Search();

    /**
     * @returns the location of inferred cghangepoints in the data . 
     */
    std::vector<int> Retreive_changepoints();

    /**
     * @returns a vector that contains the total number of intervals forming the area of life of the different candidates before the pruning step at each iteration.
     */
    std::vector<int> Get_intervals();


    /**
     * @returns a vector that contains the total number of candidates before the pruning step at each iteration.
     */
    std::vector<int> Get_candidates();

    std::vector<double> Retreive_means();
    std::vector<double> Retreive_costs();
    
};
#endif
