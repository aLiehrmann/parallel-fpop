#ifndef DEF_CANDIDATE
#define DEF_CANDIDATE

#include "Interval.h"
#include <vector>
#include "Quadratic.h"
#include <list>
#include "Ordered_list_of_intervals.h"

class Candidate
{
private:
    double cost_up_to_tau;
    int tau;
    Ordered_list_of_intervals z;
    Quadratic quad;
    double pen;

public:
    /**
     * @details Instantiates a candidate.
     * @param[in] tau_ position of the candidate's the last changepoint.
     * @param[in] z_ the sorted list of intervals forming the candidate's area of life.
     * @param[in] cost_up_to_tau_ cost of the best segmentation of the data until the point tau.
     * @param[in] pen_ penalty which depending of the length of the segment formed by the points located after tau.
     * @param[in] quad_ quadratic form build from the point located after tau.
     */
    Candidate(int tau_, Ordered_list_of_intervals z_, double cost_up_to_tau_, double pen_, Quadratic quad_);
    

    /**
     * @returns the minimum of the cost function.
     */
    double Minimum_of_cost_function();
    
    /**
     * @details updates the penalty which depending of the length of the segment formed by the points located after tau.
     * @param[in] pen_ new penalty.
     */
    void Set_penalty(double pen_);

    /**
     * @details Searches for the intervals over which the cost function of the last candidate is beaten by the cost function of past candidates.
     * Updates the area of life of the cost function of the last candidate with the complementary of the union of the intervals found.
     * @param[in] vector_of_it_candidates vector of iterator that points on candidates.
     * @param[in] D default of the last candidate's area of life.
     */
    void Compare_to_past_candidates (std::vector<std::list<Candidate>::iterator> & vector_of_it_candidates, Interval & D);

    /**
     * @details Searches the intervals over which the cost function of the current candidate is beaten by the cost function of future candidates.
     * Updates the current candidate's area of life by intersecting it with the intersection of the intervals found.
     * @param[in] vector_of_it_candidates vector of iterator that points on candidates.
     * @param[in] chosen_candidates indices of the candidates chosen for the comparison.
     */
    void Compare_to_future_candidates (std::vector<std::list<Candidate>::iterator> & vector_of_it_candidates, std::vector<int> & chosen_candidates, Interval & D);
    
    /**
     * @details Updates the quadratic form of the current candidate by adding to it the wt-weighted quadratic (y-mu)^2.
     * @param[in] wt a weight.
     * @param[in] y a point.
     */
    void Add_quadratic(double wt, double y);
    
    /**
     * @returns the sorted list of intervals forming the current candidate's area of life.
     */
    Ordered_list_of_intervals GetZ();

    /**
     * @returns tau of the current candidate.
     */
    int Get_tau();

    double Argmin_of_cost_function();
    void Compare_to_last_candidates (Candidate & last_candidate, Interval & D);
};


#endif