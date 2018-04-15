/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef UNITS_HPP_INCLUDED
#define UNITS_HPP_INCLUDED

#pragma GCC system_header

#ifdef NO_UNITS
#  include "impl/no_units.hpp"
#else
#  include "impl/with_units.hpp"
#endif // NO_UNITS

namespace cellsim {
     using units_impl::concentration_ut;
     using units_impl::electric_potential_ut;
     using units_impl::potential_SR_ut;
     using units_impl::length_ut;
     using units_impl::flux_ut;
     using units_impl::time_ut;
     using units_impl::inverse_time_ut;
     using units_impl::diff_coef_ut;
     using units_impl::permeablity_ut;

     using units_impl::conductance_ut;
     using units_impl::gamma_ut;
     using units_impl::beta_ut;
} // namespace cellsim

// clean up
#undef MAKE_LITERAL_OPERATOR
#undef MAKE_LITERAL_OPERATOR_BOOST

#endif //UNITS_HPP_INCLUDED
