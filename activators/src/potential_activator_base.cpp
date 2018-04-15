#include "potential_activator_base.hpp"

cellsim::Potential_Activator_Base::Potential_Activator_Base(
     electric_potential_ut v_act,
     time_ut t0,
     double alpha)
     : v0_(0_mV),
       v_act_(v_act),
       alpha_(alpha),
       t_(0.),
       t0_(t0.value())
{}
