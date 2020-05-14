#include <vector>
#include <iostream>
#include "Fpop.h"
#include "Candidate.h"
#include <algorithm>
#include <math.h>
#include <list>
#include <limits>
#include <atomic>
#include "Ordered_list_of_intervals.h"
#include "Interval.h"

#include "omp.h"


//####### Constructor #######////####### Constructor #######////####### Constructor #######//
//####### Constructor #######////####### Constructor #######////####### Constructor #######//


Fpop::Fpop(std::vector<double> y_,
    double lambda_,
    double muMinLocal,
    double muMaxLocal,
    std::vector<double>  wt_)
{
    y = y_;
    n = y_.size();
    d = Interval(muMinLocal, muMaxLocal);

    y.insert(y.begin(), 0);
    lambda = lambda_;

    if (wt_.size() == 1 && wt_[0] == 0)
    {
        wt = std::vector<double> (y.size(), 1);
    } else {
       wt = wt_;
       wt.insert(wt.begin(), 0);
    }
    cp = std::vector<int> (y.size(), 0);
    costs = std::vector<double> (y.size(), 0);
    means = std::vector<double> (y.size(), 0);
    nb_candidates = std::vector<int> (y.size()-1, 0);
    nb_intervals = std::vector<int> (y.size()-1, 0);
}

Fpop::Fpop(){}


//####### changepoints_search #######////####### changepoints_search #######////####### changepoints_search #######//
//####### changepoints_search #######////####### changepoints_search #######////####### changepoints_search #######//


void Fpop::Search()
{
    std::list<Candidate> list_of_candidates {Candidate(0,  Ordered_list_of_intervals (d), 0, 0, Quadratic())};
    double F;
    double ARG_F;
    int t_hat;
    double min_candidate;
    double arg_min_candidate;
    int index;
    std::vector<int> chosen_candidates;
    std::vector<std::list<Candidate>::iterator> vector_of_it_candidates;

    for (int t {1}; t < y.size(); t++)
    {
        /*
            We initialize a vector of iterators.
            Each iterator points to a candidate in list_of_candidates.
            At this step, the last element of the vector does not point to any candidate.
        */

        std::vector<std::list<Candidate>::iterator> vector_of_it_candidates (list_of_candidates.size() + 1);
        index = 0;
        for (auto it_candidate {list_of_candidates.begin()}; it_candidate != list_of_candidates.end(); ++it_candidate)
        {
            vector_of_it_candidates[index] = it_candidate;
            index += 1;
        }

        /*
            (1) The quadratic form is updated for all the candidates still considered.

            (2) The length-dependent penalty for the last segment is updated for all candidates still considered.

            (3) we update the minimum cost of segmentation and associated candidate.
        */

        F = std::numeric_limits<double>::max();
        for (int i {0}; i < vector_of_it_candidates.size() - 1; i++)
        {
            (*vector_of_it_candidates[i]).Add_quadratic(wt[t], y[t]); //1
            min_candidate = (*vector_of_it_candidates[i]).Minimum_of_cost_function();
            arg_min_candidate = (*vector_of_it_candidates[i]).Argmin_of_cost_function();
            if (min_candidate < F) //(3)
            {
                F = min_candidate;
                ARG_F = arg_min_candidate;
                t_hat = (*vector_of_it_candidates[i]).Get_tau();
            }
        }

        /*
            (1) We save the position of the last changepoint of the candidate that minimizes the cost of segmentation up to the point t.

            (2) A new candidate is introduced whose last changepoint corresponds to point t.

            (3) The last element of array_of_candidates is now pointing to the last introduced candidate.
        */

        cp[t] = t_hat; //(1)
        costs[t] = F;
        means[t] = ARG_F;
        list_of_candidates.push_back( Candidate(t, Ordered_list_of_intervals (d), F + lambda, 0, Quadratic())); //(2)
        vector_of_it_candidates[vector_of_it_candidates.size() - 1] = --list_of_candidates.end(); //(3)

        /*
           (1) We save the sum of the intervals in candidates's area of life.

           (2) We save the number of candidates still considered.
        */

        for (int i {0}; i < vector_of_it_candidates.size(); i++)
        {
            nb_intervals[t - 1] += (*vector_of_it_candidates[i]).GetZ().size(); //(1)
        }
        nb_candidates[t - 1] += list_of_candidates.size(); //(2)

        // We update the the last candidate's area of life.

        (*vector_of_it_candidates.back()).Compare_to_past_candidates(vector_of_it_candidates, d);

        // Candidates whose area of life is empty are pruned.

        list_of_candidates.erase(std::remove_if(list_of_candidates.begin(), list_of_candidates.end(), [](Candidate & a) {
            return a.GetZ().Is_empty();
        }), list_of_candidates.end());
    }
}


void Fpop::Search_parallel(int tid, int nbThreads, double * F, double * ARG_F, int * t_hat,
    bool * firstMin, double * F_min, double * ARG_F_min, int * t_hat_min)
{
    std::list<Candidate> list_of_candidates {Candidate(0,  Ordered_list_of_intervals (d), 0, 0, Quadratic())};
    // double F;
    // int t_hat;
    double min_candidate;
    double arg_min_candidate;
    // double ARG_F;
    int index;
    std::vector<int> chosen_candidates;
    std::vector<std::list<Candidate>::iterator> vector_of_it_candidates;


    for (int t {1}; t < y.size(); t++)
    {
        /*
            We initialize a vector of iterators.
            Each iterator points to a candidate in list_of_candidates.
            At this step, the last element of the vector does not point to any candidate.
        */

        std::vector<std::list<Candidate>::iterator> vector_of_it_candidates (list_of_candidates.size() + 1);
        index = 0;
        for (auto it_candidate {list_of_candidates.begin()}; it_candidate != list_of_candidates.end(); ++it_candidate)
        {
            vector_of_it_candidates[index] = it_candidate;
            index += 1;
        }

        /*
            (1) The quadratic form is updated for all the candidates still considered.
            
            (2) The length-dependent penalty for the last segment is updated for all candidates still considered.
            
            (3) we update the minimum cost of segmentation and associated candidate.
        */

        F[tid] = std::numeric_limits<double>::max();
        for (int i {0}; i<vector_of_it_candidates.size() - 1; i++)
        {
            (*vector_of_it_candidates[i]).Add_quadratic(wt[t], y[t]); //1
            min_candidate = (*vector_of_it_candidates[i]).Minimum_of_cost_function();
            arg_min_candidate = (*vector_of_it_candidates[i]).Argmin_of_cost_function();
            if (min_candidate < F[tid]) //(3)
            {
                F[tid] = min_candidate;
                ARG_F[tid] = arg_min_candidate;
                t_hat[tid] = (*vector_of_it_candidates[i]).Get_tau();
            }
        }

        //TODO Reduce F between the threads, and give F, ARG_F, and t_hat corresponding to the smallest F
        #pragma omp critical
        {
            if (!(*firstMin))
            {
                (*F_min) = F[tid];
                (*ARG_F_min) = ARG_F[tid];
                (*t_hat_min) = t_hat[tid];
                (*firstMin) = true;
            } else {
                if (F[tid] < (*F_min))
                {
                    (*F_min) = F[tid];
                    (*ARG_F_min) = ARG_F[tid];
                    (*t_hat_min) = t_hat[tid];
                }
            }
        }
        #pragma omp barrier

        #pragma omp single
        {
            (*firstMin) = false;
        } // Implicit barrier

        // #pragma omp barrier
        // #pragma omp single
        // {
        //     double min_F = F[0];
        //     int min_tid = 0;
        //
        //     for (int i = 1; i < nbThreads; ++i)
        //     {
        //         if (F[i] < min_F)
        //         {
        //             min_F = F[i];
        //             min_tid = i;
        //         }
        //     }
        //
        //     for (int i = 0; i < nbThreads; ++i)
        //     {
        //         F[i] = min_F;
        //         ARG_F[i] = ARG_F[min_tid];
        //         t_hat[i] = t_hat[min_tid];
        //     }
        // } // Implicit barrier

        /*
            (1) We save the position of the last changepoint of the candidate that minimizes the cost of segmentation up to the point t.
            
            (2) A new candidate is introduced whose last changepoint corresponds to point t.
            
            (3) The last element of array_of_candidates is now pointing to the last introduced candidate.
        */

        // cp[t] = t_hat[tid]; //(1)
        // costs[t] = F[tid];
        // means[t] = ARG_F[tid];
        cp[t] = (*t_hat_min); //(1)
        costs[t] = (*F_min);
        means[t] = (*ARG_F_min);
        // list_of_candidates.push_back( Candidate(t, Ordered_list_of_intervals (d), F[tid]+lambda, 0, Quadratic())); //(2)
        list_of_candidates.push_back( Candidate(t, Ordered_list_of_intervals (d), (*F_min) + lambda, 0, Quadratic())); //(2)
        vector_of_it_candidates[vector_of_it_candidates.size()-1] = --list_of_candidates.end(); //(3)

        /*
           (1) We save the sum of the intervals in candidates's area of life.
            
           (2) We save the number of candidates still considered.
        */

        for (int i {0}; i < vector_of_it_candidates.size(); i++)
        {
            nb_intervals[t - 1] += (*vector_of_it_candidates[i]).GetZ().size(); //(1)
        }
        nb_candidates[t - 1] += list_of_candidates.size(); //(2)

        // We update the area of life of the candidates.

        (*vector_of_it_candidates.back()).Compare_to_past_candidates(vector_of_it_candidates, d);

        //Candidates whose area of life is empty are pruned.

        list_of_candidates.erase(std::remove_if(list_of_candidates.begin(), list_of_candidates.end(), [](Candidate & a) {
            return a.GetZ().Is_empty();
        }), list_of_candidates.end());

     }

}


void Fpop::Search_parallel_loops(int nbThreads)
{
    std::list<Candidate> list_of_candidates {Candidate(0,  Ordered_list_of_intervals (d), 0, 0, Quadratic())};
    double F;
    double ARG_F;
    int t_hat;
    double min_candidate;
    double arg_min_candidate;
    int index;
    std::vector<int> chosen_candidates;
    std::vector<std::list<Candidate>::iterator> vector_of_it_candidates;

    std::atomic<bool> computationDone(false);
    std::atomic<int> nbJobs(0);
    std::atomic<int> nbThreadsDone(0);

    std::list<Interval> list_of_intervals_array[nbThreads];
    int globalLoopSize;

    #pragma omp parallel num_threads(nbThreads)
    {
        int tid = omp_get_thread_num();

        if (0 == tid)
        {
            // Main computing thread
            for (int t {1}; t < y.size(); t++)
            {
                /*
                    We initialize a vector of iterators.
                    Each iterator points to a candidate in list_of_candidates.
                    At this step, the last element of the vector does not point to any candidate.
                */

                std::vector<std::list<Candidate>::iterator> vector_of_it_candidates (list_of_candidates.size() + 1);
                index = 0;
                for (auto it_candidate {list_of_candidates.begin()}; it_candidate != list_of_candidates.end(); ++it_candidate)
                {
                    vector_of_it_candidates[index] = it_candidate;
                    index += 1;
                }

                /*
                    (1) The quadratic form is updated for all the candidates still considered.

                    (2) The length-dependent penalty for the last segment is updated for all candidates still considered.

                    (3) we update the minimum cost of segmentation and associated candidate.
                */

                F = std::numeric_limits<double>::max();
                for (int i {0}; i < vector_of_it_candidates.size() - 1; i++)
                {
                    (*vector_of_it_candidates[i]).Add_quadratic(wt[t], y[t]); //1
                    min_candidate = (*vector_of_it_candidates[i]).Minimum_of_cost_function();
                    arg_min_candidate = (*vector_of_it_candidates[i]).Argmin_of_cost_function();
                    if (min_candidate < F) //(3)
                    {
                        F = min_candidate;
                        ARG_F = arg_min_candidate;
                        t_hat = (*vector_of_it_candidates[i]).Get_tau();
                    }
                }

                /*
                    (1) We save the position of the last changepoint of the candidate that minimizes the cost of segmentation up to the point t.

                    (2) A new candidate is introduced whose last changepoint corresponds to point t.

                    (3) The last element of array_of_candidates is now pointing to the last introduced candidate.
                */

                cp[t] = t_hat; //(1)
                costs[t] = F;
                means[t] = ARG_F;
                list_of_candidates.push_back( Candidate(t, Ordered_list_of_intervals (d), F + lambda, 0, Quadratic())); //(2)
                vector_of_it_candidates[vector_of_it_candidates.size() - 1] = --list_of_candidates.end(); //(3)

                /*
                   (1) We save the sum of the intervals in candidates's area of life.

                   (2) We save the number of candidates still considered.
                */

                for (int i {0}; i < vector_of_it_candidates.size(); i++)
                {
                    nb_intervals[t - 1] += (*vector_of_it_candidates[i]).GetZ().size(); //(1)
                }
                nb_candidates[t - 1] += list_of_candidates.size(); //(2)

                // We update the the last candidate's area of life.

                globalLoopSize = vector_of_it_candidates.size();
                nbJobs++;

                int sizePerThread = globalLoopSize / nbThreads;
                int threadsWithOneMore = globalLoopSize - (sizePerThread * nbThreads);
                int beginLoop;
                int endLoop;

                // Ensure that we are not going out of bounds
                if (tid < threadsWithOneMore)
                {
                    beginLoop = tid * sizePerThread + tid;
                    endLoop = beginLoop + sizePerThread + 1;
                } else {
                    beginLoop = tid * sizePerThread + threadsWithOneMore;
                    endLoop = beginLoop + sizePerThread;
                }

                // Compute a fraction of the loop
                list_of_intervals_array[tid] = (*vector_of_it_candidates.back()).Compare_to_past_candidates_parallel(vector_of_it_candidates, d, beginLoop, endLoop);
                // (*vector_of_it_candidates.back()).Compare_to_past_candidates(vector_of_it_candidates, d);

                // The main thread is done with its part of the computation
                nbThreadsDone++;

                while (nbThreadsDone < nbThreads)
                {
                    // Wait for all the threads to finish their part of the work
                }

                Ordered_list_of_intervals list_of_merged_intervals (list_of_intervals_array[tid]);
                for (int i = 1; i < nbThreads; ++i)
                {
                    list_of_merged_intervals.Append(list_of_intervals_array[i]);
                }
                list_of_merged_intervals.Complementary_in(d);
                (*vector_of_it_candidates.back()).SetZ(list_of_merged_intervals);

                // Candidates whose area of life is empty are pruned.

                list_of_candidates.erase(std::remove_if(list_of_candidates.begin(), list_of_candidates.end(), [](Candidate & a) {
                    return a.GetZ().Is_empty();
                }), list_of_candidates.end());
            }

            // The total computation is done, so indicate it to the other threads
            computationDone = true;
        } else {
            // Threads that "help" to compute some of the loops
            int nbJobsLocal = 0;
            while (!computationDone)
            {
                while (nbJobsLocal == nbJobs)
                {
                    // While loop to wait for the main thread to attribute some work
                }

                nbJobsLocal++;

                int sizePerThread = globalLoopSize / nbThreads;
                int threadsWithOneMore = globalLoopSize - (sizePerThread * nbThreads);
                int beginLoop;
                int endLoop;

                // Ensure that we are not going out of bounds
                if (tid < threadsWithOneMore)
                {
                    beginLoop = tid * sizePerThread + tid;
                    endLoop = beginLoop + sizePerThread + 1;
                } else {
                    beginLoop = tid * sizePerThread + threadsWithOneMore;
                    endLoop = beginLoop + sizePerThread;
                }

                // Compute a fraction of the loop
                list_of_intervals_array[tid] = (*vector_of_it_candidates.back()).Compare_to_past_candidates_parallel(vector_of_it_candidates, d, beginLoop, endLoop);

                // Indicate to the main thread that this current thread is done for this computation
                nbThreadsDone++;
            }
        }
    }
}


//####### Retreive_changepoints #######////####### Retreive_changepoints #######////####### Retreive_changepoints #######//
//####### Retreive_changepoints #######////####### Retreive_changepoints #######////####### Retreive_changepoints #######//


std::vector<int> Fpop::Retrieve_changepoints()
{
    std::vector<int> list_of_changepoints;
    int i (y.size()-1);
    while (cp[i] != 0)
    {
        list_of_changepoints.push_back(cp[i]);
        i = cp[i];
    }
    std::reverse(list_of_changepoints.begin(), list_of_changepoints.end());
    list_of_changepoints.push_back(y.size()-1);
    return list_of_changepoints;
}

//####### Retreive_means #######////####### Retreive_means #######////####### Retreive_means #######//
//####### Retreive_means #######////####### Retreive_means #######////####### Retreive_means #######//


std::vector<double> Fpop::Retrieve_means()
{
    std::vector<double> list_of_means;
    int i (y.size()-1);
    while (cp[i] != 0)
    {
        list_of_means.push_back(means[i]);
        i = cp[i];
    }
    list_of_means.push_back(means[i-1]);
    std::reverse(list_of_means.begin(), list_of_means.end());
    return list_of_means;
}

//####### Retreive_costs #######////####### Retreive_costs #######////####### Retreive_costs #######//
//####### Retreive_costs #######////####### Retreive_costs #######////####### Retreive_costs #######//

std::vector<double> Fpop::Retrieve_costs()
{
    std::vector<double> list_of_costs;
    int i (y.size()-1);
    while (cp[i] != 0)
    {
        list_of_costs.push_back(costs[i]);
        i = cp[i];
    }
    list_of_costs.push_back(costs[i-1]);
    std::reverse(list_of_costs.begin(), list_of_costs.end());
    return list_of_costs;
}


//####### number_of_candidates_in_list #######////####### number_of_candidates_in_list #######////####### number_of_candidates_in_list #######//
//####### number_of_candidates_in_list #######////####### number_of_candidates_in_list #######////####### number_of_candidates_in_list #######//


std::vector<int> Fpop::Get_candidates()
{
    return nb_candidates;
}


//####### number_of_intervals_in_all_areas_of_life #######////####### number_of_intervals_in_all_areas_of_life #######////####### number_of_intervals_in_all_areas_of_life #######//
//####### number_of_intervals_in_all_areas_of_life #######////####### number_of_intervals_in_all_areas_of_life #######////####### number_of_intervals_in_all_areas_of_life #######//


std::vector<int> Fpop::Get_intervals()
{
    return nb_intervals;
}
