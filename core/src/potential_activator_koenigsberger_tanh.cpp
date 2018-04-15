#include "potential_activator_koenigsberger_tanh.hpp"

#include "cell_koenigsberger.hpp"

#include <cmath>

namespace activator = cellsim::model_koenigsberger::activator;

activator::Potential_Activator_Tanh::Potential_Activator_Tanh(
     size_type start_column,
     size_type end_column,
     electric_potential_ut v_act,
     time_ut t0,
     double alpha)
     : Potential_Activator(start_column, end_column, v_act, t0, alpha),
       t_(-4. / alpha), 
       t_fin_(4. / alpha)
{}

void activator::Potential_Activator_Tanh::operator() (Cell& c, time_ut dt)
{
     if (v0_ == 0_mV) {
	  v0_ = c.v();
     }
     
     if (t_ < t_fin_) {
	  c.v_new() = v0_ + 0.5 * v_act_ * (std::tanh(alpha_ * t_)+1);
     }
     else {
	  // After Tanh increase, holding v_act_ for t0_
	  // (v_act_ is modified just below for final exponential decay)
	  c.v_new() = v0_ + v_act_;	  
     }
     
     // Exponential decay after holding
     if ((t_fin_ + t0_) < t_) {
	  v_act_ *= (1. - alpha_ * dt.value());
     }
     t_ += dt.value();
}

