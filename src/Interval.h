#ifndef DEF_INTERVAL
#define DEF_INTERVAL

#include <vector>
#include <array> 
#include <list>

class Interval
{
private:
    double begin;
    double end;
public:
    /**
     * @details Instanciates an empty interval. 
     * Any interval of the form [a, b] such that a>=b is an empty interval (by default [+ 1, -1]).
     */
    Interval();

    /**
     * @details Instanciates an interval of the form [begin_,end_].
     */
    Interval(double begin_, double end_);

    /**
     * @details Instantiates the intersection of the intervals contained in a list of intervals. This intersection may be empty or a singleton.
     * @param[in] list_of_intervals_to_intersect an unsorted list of nonempty intervals, which are not singletons. This list can be empty, in which case it instantiates an empty interval.
     */
    Interval(std::list<Interval> & list_of_intervals_to_intersect);

    /**
     * @details Updates the current interval by intersecting it with another interval.
     * @param[in] interval_to_intesect an interval that can be empty or take the form of a singleton.
     */
    void operator&=(Interval const& interval_to_intesect);

    /**
     * @details Compare the beginning of two intervals.
     * @returns TRUE si le début de l'intervalle 'interval1' est inférieur au début de l'intervalle 'interval2', FALSE otherwise.
     */
    static bool Compare_begin(Interval const& interval1, Interval const& interval2);

    /**
     * @details Compare the end of the two intervals.
     * @returns TRUE if the beginning of the interval 'interval1' is less than the beginning of the interval 'interval2', FALSE otherwise.
     */
    static bool Compare_end(Interval const& interval1, Interval const& interval2);

    /**
     * @returns TRUE if the current interval is empty ([a, b] such that a> b) or is a singleton ([a, a]), FALSE otherwise. 
     */
    bool IsEmpty_or_singleton();

    /**
     * @returns the beginning of the current interval.
     */
    double Get_begin();

    /**
     * @returns the end of the current interval.
     */
    double Get_end();

    /**
     * @details updates the beginning of the current interval.
     */
    void Set_begin(double begin_);

    /**
     * @details updates the end of the current interval.
     */
    void Set_end(double end_);

    
};

#endif