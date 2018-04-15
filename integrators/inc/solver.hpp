/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED

#include "definitions.hpp"

namespace cellsim {

     GCC_DIAG_OFF(float-equal)

     template <typename vector_type>
     void triangular_solve(const vector_type a, const vector_type b, const vector_type c,
			   const vector_type rhs, vector_type& u)
     //Solves for a vector u[1..n] the tridiagonal linear set given by equation (2.4.1). a[1..n], b[1..n], c[1..n], and r[1..n] are input vectors and are not modified.
     {
	  const size_type N(b.size());
	  
	  // Check dimensions are compatible
	  assert(!a.empty() && a.size() == c.size() && a.size() == (b.size()-1) &&
		 u.size() == b.size() && rhs.size() == b.size());

	  typename vector_type::value_type bet(b[1]);
	  vector_type gam;

	  gam = vector_type(b.size());
	  if (bet == 0.0) {
	       std::cerr << "Error 1 in tridag\n";
	       return;
	  }
	  // If this happens then you should rewrite your equations as a set of order N − 1, with u2 trivially eliminated.
     
	  u[1] = rhs[1] / bet;

	  for (size_type j(2) ; j <= N ; ++j) {     //Decompositiona nd forward substitution.
	       gam[j] = c[j-1] / bet;
	       bet = b[j] - a[j]*gam[j];

	       if (bet == 0.0) {
		    std::cerr << "Error 2 in tridag\n";
		    return;
	       }

	       u[j] = (rhs[j] - a[j]*u[j-1]) / bet;
	  }

	  for (size_type j(N-1) ; j >= 1 ; --j) {
	       u[j] -= gam[j+1] * u[j+1]; 
	  }
     }

     GCC_DIAG_ON(float-equal)
}

#endif //SOLVER_HPP_INCLUDED
