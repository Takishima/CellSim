/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef RICHARDSON_EXTRAPOLATOR_HPP_INCLUDED
#define RICHARDSON_EXTRAPOLATOR_HPP_INCLUDED

#include "equation_type.hpp"

template <typename integrator_t>
class Richardson_Extrapolator
{
public:
     template <typename number_t,
	       typename variable_t>
     static void extrapolate(const number_t y_now, 
			     number_t& y_new
			     const variable_t a,
			     const variable_t da)
	  {
	  }
     
private:
};

#endif //RICHARDSON_EXTRAPOLATOR_HPP_INCLUDED
