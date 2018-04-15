/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef FTCS_IMPL_HPP_INCLUDED
#define FTCS_IMPL_HPP_INCLUDED

#include <cassert>
#include <functional>
#include <utility>
#include <vector>
#include <iostream>

#include "config.hpp"
#include "equation_type.hpp"
#include "no_op.hpp"
#include "point_container.hpp"
#include "units.hpp"
#include "apply_activator.hpp"

namespace cellsim {
namespace integrators {
namespace ftcs_impl_ {
     struct default_t;

     template <typename T>
     struct type_helper
     {
	  typedef time_ut time_t;
	  typedef length_ut length_t;
	  typedef diff_coef_ut diff_coef_t;
     };
     
     template <>
     struct type_helper<double>
     {
	  typedef double time_t;
	  typedef double length_t;
	  typedef double diff_coef_t;
     };

} // namespace ftcs_impl_
} // namespace integrators
} // namespace cellsim

#endif //FTCS_IMPL_HPP_INCLUDED
