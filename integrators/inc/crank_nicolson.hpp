/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef CRANK_NICOLSON_HPP_INCLUDED
#define CRANK_NICOLSON_HPP_INCLUDED

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "equation_type.hpp"
#include "solver.hpp"

namespace cellsim {
namespace integrators {

     /*!
      * \brief Integrator for the diffusion equation using the Crank-Nicolson scheme.
      *
      * The equation we want to integrate is of the following form : \f[\frac{\partial f}{\partial t} = D \frac{\partial^2 f}{\partial x^2} + k f\f]
      * where \f$ k \f$ is a constant.
      */
     class CN_Integrator
     {
     public:
	  typedef std::vector<double> vector_type;
	  typedef vector_type::value_type value_type;
	  typedef vector_type::size_type size_type;
	  
	  /*!
	   * Simple constructor
	   *
	   * \param N Size of system
	   * \param D Diffusion coefficient
	   * \param k REMOVE ?
	   * \param dt Time step
	   * \param dx Spatial step
	   */
	  CN_Integrator(size_type N, double D, double k,
			double dt, double dx);
	  
	  void setup_system(double D, double k, double dt, double dx);
	  void setup_system(size_type N, double D, double k,
			    double dt, double dx);

	  template <typename read_func, typename write_func, typename array_type>
	  vector_type do_step(const double t,
			      array_type& y_now,
			      const read_func read_mem_func,
			      const write_func write_mem_func) const;
     private:
	  double r_;
	  double k_;

	  vector_type a_;
	  vector_type b_;
	  vector_type c_;
     };

     /*!
      * Integrate the diffusion equation one time step.
      * 
      *
      */
     template <typename read_func, typename write_func, typename array_type>
     CN_Integrator::vector_type 
     CN_Integrator::do_step(const double /*t*/,
			    array_type& y_now,
			    const read_func read_mem_func,
			    const write_func write_mem_func) const
     {
	  assert(y_now.size() == b_.size());

	  vector_type u(y_now.size(), value_type());
	  vector_type rhs;
	  rhs.reserve(y_now.size());

	  // Build RHS of system to solve
	  const auto END(end(y_now));
	  for (auto prev(END), it(begin(y_now)), next(begin(y_now)+1) ;
	       it != END ; ++it, ++next) {

	       value_type tmp(0);
	       
	       if (prev != END) {
		    //std::cout << "p = " << read_mem_func(*prev) << " ";
		    tmp += r_ * (read_mem_func(*prev));
		    ++prev;
	       }
	       else {
		    prev = begin(y_now);
	       }

	       //std::cout << "it = " << read_mem_func(*it) << " ";
	       tmp += (1 - 2*r_) * (read_mem_func(*it));

	       if (next != END) {
		    //std::cout << "n = " << read_mem_func(*next) << " ";
		    tmp += r_ * (read_mem_func(*next));
	       }
	       //std::cout << std::endl;
	       rhs.push_back(tmp);
	  }
	  
	  triangular_solve(a_, b_, c_, rhs, u);

	  const auto END2 (end(u));
	  auto y_it(begin(y_now));

	  // Store update original array with values at new timestep
	  for (auto u_it (begin(u)) ; u_it != END2 ; ++u_it, ++y_it) {
	       write_mem_func(*y_it, *u_it);
	  }
	  
	  return u;
     }
	 
} // namespace integrators
} // namespace cellsim

#endif //CRANK_NICOLSON_HPP_INCLUDED
