#ifndef DEF_ORDERED_LIST_OF_INTERVALS
#define DEF_ORDERED_LIST_OF_INTERVALS

#include <list>
#include "Interval.h"


class Ordered_list_of_intervals
{
private:
    std::list<Interval> list_of_intervals;
public:
    /**
     * @details Instantiates the sorted union of the intervals contained in a list of intervals.
     * @param[in] list_of_intervals_to_merge a list of nonempty intervals, not containing singletons, unordered. This list can be empty.
     */
    Ordered_list_of_intervals(std::list<Interval> & list_of_intervals_to_merge);

    void Append(std::list<Interval> & list_of_intervals_to_merge);

    /**
     * @details Instantiates an ordered list of intervals with interval 'interval'.
     */
    Ordered_list_of_intervals(Interval interval);

    /**
     * @details Instantiates an empty ordered list of intervals.
     */
    Ordered_list_of_intervals();

    /**
     * @details Updates the current ordered list of intervals by intersecting each interval contained in this list with an interval 'interval_to_intersect'.
     */
    void Intersect_with(Interval & interval_to_intersect);

    /**
     * @details Updates the ordered list of intervals with its complementary in the interval 'D'. If the list is empty, it is updated with the interval D.
     */
    void Complementary_in(Interval & d);

    /**
     * @returns TRUE if the ordered list of intervals is empty, false otherwise.
     */
    bool Is_empty();

    /**
     * @returns the length of the ordered list of intervals.
     */
    int size();

    /**
     * @returns the ordered list of intervals.
     */
    std::list<Interval> Get_list();
};


#endif
