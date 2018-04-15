#include "point_index.hpp"

cellsim::Point_Index::Point_Index(const Point_Container* container, 
				  size_type index_a)
     : index_(index_a + 1), container_(container)
{
     /*
      * Note: we must initialize index_ at index + 1
      *       (see Point_Container's constructor for more info)
      */
}
