#include "Candidate.h"
#include "Quadratic.h"
#include "Interval.h"
#include <vector>
#include <array>
#include <iostream>
#include <list>
#include <math.h>
#include "Ordered_list_of_intervals.h"


//####### Constructor #######////####### Constructor #######////####### Constructor #######//
//####### Constructor #######////####### Constructor #######////####### Constructor #######//


Candidate::Candidate(int tau_, Ordered_list_of_intervals z_, double cost_up_to_tau_, double pen_, Quadratic quad_)
{

    tau = tau_;
    z = z_;
    cost_up_to_tau = cost_up_to_tau_;
    pen = pen_;
    quad = quad_;
}


//####### Minimum_of_cost_function #######////####### Minimum_of_cost_function #######////####### Minimum_of_cost_function #######//
//####### Minimum_of_cost_function #######////####### Minimum_of_cost_function #######////####### Minimum_of_cost_function #######//


double Candidate::Minimum_of_cost_function()
{
    return quad.Minimum() + pen + cost_up_to_tau;
}

double Candidate::Argmin_of_cost_function()
{
    return quad.Argmin();
}


//####### update_penalty_on_last_segment #######////####### update_penalty_on_last_segment #######////####### update_penalty_on_last_segment #######//
//####### update_penalty_on_last_segment #######////####### update_penalty_on_last_segment #######////####### update_penalty_on_last_segment #######//


void Candidate::Set_penalty(double pen_)
{
    pen = pen_;
}


//####### add_new_point #######////####### add_new_point #######////####### add_new_point #######//
//####### add_new_point #######////####### add_new_point #######////####### add_new_point #######//


void Candidate::Add_quadratic(double wt, double y)
{
    quad.Add_coef(wt*pow(y,2), -2*wt*y, wt);
}


//####### Compare_to_past_candidate #######////####### Compare_to_past_candidate #######////####### Compare_to_past_candidate #######//
//####### Compare_to_past_candidate #######////####### Compare_to_past_candidate #######////####### Compare_to_past_candidate #######//


void Candidate::Compare_to_past_candidates (std::vector<std::list<Candidate>::iterator> & vector_of_it_candidates, Interval & D)
{
    std::list<Interval> list_of_intervals;
    Interval interval;
    Quadratic new_quad;
    for (int i {0}; i < vector_of_it_candidates.size() - 1; i++)
    {

        new_quad = (*vector_of_it_candidates[i]).quad - quad;
        new_quad.Add_coef((*vector_of_it_candidates[i]).cost_up_to_tau - cost_up_to_tau, 0, 0);
        interval = new_quad.Negative_interval(D);
        if (!interval.IsEmpty_or_singleton())
        {
            list_of_intervals.push_back(interval);
        }
        vector_of_it_candidates[i]->z.Intersect_with(interval);

    }
    Ordered_list_of_intervals list_of_merged_intervals (list_of_intervals);
    list_of_merged_intervals.Complementary_in(D);
    z = list_of_merged_intervals;
}


std::list<Interval> Candidate::Compare_to_past_candidates_parallel (std::vector<std::list<Candidate>::iterator> & vector_of_it_candidates, Interval & D, int begin, int end)
{
    std::list<Interval> list_of_intervals;
    Interval interval;
    Quadratic new_quad;
    for (int i = begin; i < end; i++)
    {
        new_quad = (*vector_of_it_candidates[i]).quad - quad;
        new_quad.Add_coef((*vector_of_it_candidates[i]).cost_up_to_tau - cost_up_to_tau, 0, 0);
        interval = new_quad.Negative_interval(D);
        if (!interval.IsEmpty_or_singleton())
        {
            list_of_intervals.push_back(interval);
        }
        vector_of_it_candidates[i]->z.Intersect_with(interval);
    }
    // Ordered_list_of_intervals list_of_merged_intervals (list_of_intervals);
    // list_of_merged_intervals.Complementary_in(D);
    // z = list_of_merged_intervals;
    return list_of_intervals;
}


//####### Get_last_changepoint_position #######////####### Get_last_changepoint_position #######////####### Get_last_changepoint_position #######//
//####### Get_last_changepoint_position #######////####### Get_last_changepoint_position #######////####### Get_last_changepoint_position #######//

int Candidate::Get_tau()
{
    return tau;
}


//####### get_candidate_area_of_life #######////####### get_candidate_area_of_life #######////####### get_candidate_area_of_life #######//
//####### get_candidate_area_of_life #######////####### get_candidate_area_of_life #######////####### get_candidate_area_of_life #######//


Ordered_list_of_intervals Candidate::GetZ()
{
    return z;
}


void Candidate::SetZ(Ordered_list_of_intervals & list_of_intervals)
{
    z = list_of_intervals;
}
