#include "PLC_activator.hpp"

namespace activator = cellsim::model_halidi::activator;

activator::PLC_Activator::PLC_Activator(size_type cell_idx,
					size_type center_idx,
					size_type region_size,
					flux_ut J_act)
     : Activator_Base(cell_idx),
       J_act_(J_act), center_(center_idx), region_size2_(region_size/2)
{}

// cellsim::flux_ut 
// activator::PLC_Activator::operator()(const Point_Index i) const
// {
//      size_type d(0);
//      if (i.index() > center_) {
// 	  d = i.index() - center_;
//      }
//      else if (i.index() < center_) {
// 	  d = center_ - i.index() + 1;
//      }

//      if (d > region_size2_) {
// 	  return flux_ut();
//      }
//      else {
// 	  return J_act_;
//      }
// }

