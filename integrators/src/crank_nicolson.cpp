#include "crank_nicolson.hpp"

namespace integrators = cellsim::integrators;

integrators::CN_Integrator::CN_Integrator(size_type N, double D, double k,
					  double dt, double dx)
     : r_(0.), k_(0), a_(N-1), b_(N), c_(N-1)
{
     assert(N > 1);
     setup_system(D, k, dt, dx);
}

void integrators::CN_Integrator::setup_system(double D, double k, double dt, double dx)
{
     k_ = k;
     r_ = (D * dt / (2*dx*dx));

     for (auto& i : b_) {
	  i = 2*r_ + 1;
     }

     for (auto& i : a_) {
	  i = -r_;
     }

     for (auto& i : c_) {
	  i = -r_;
     }
}

void integrators::CN_Integrator::setup_system(size_type N, 
					      double D,
					      double k,
					      double dt, 
					      double dx)
{
     assert(N > 1);

     a_.resize(N-1);
     b_.resize(N);
     c_.resize(N-1);
     setup_system(D, k, dt, dx);
}
