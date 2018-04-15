#include "potassium_chloride_activator_koenigsberger.hpp"

namespace activator = cellsim::model_koenigsberger::activator;

activator::K_Cl_Activator::K_Cl_Activator(size_type start_column,
					  size_type end_column,
					  time_ut t0,
					  double alpha_rise,
					  double alpha_decay,
					  const Constants_Values& constants,
					  electric_potential_ut vK,
					  electric_potential_ut vCl)
     : Activator_Base(start_column, end_column),
       K_Cl_Activator_Base(constants, t0, alpha_rise, alpha_decay, vK, vCl)
{}
