#include "potassium_chloride_activator.hpp"

namespace activator = cellsim::model_halidi::activator;

activator::K_Cl_Activator::K_Cl_Activator(size_type cell_idx,
					  time_ut t0,
					  double alpha_rise,
					  double alpha_decay,
					  const Constants_Values& constants,
					  electric_potential_ut vK,
					  electric_potential_ut vCl)
     : Activator_Base(cell_idx),
       K_Cl_Activator_Base(constants, t0, alpha_rise, alpha_decay, vK, vCl)
{}
