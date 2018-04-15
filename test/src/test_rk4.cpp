#include <iostream>
#include <iomanip>
#include <cmath>

#include "rk4.hpp"

namespace cellsim {
namespace test {

     double f(double y, double x)
     {
	  return (-2.0 * x*y);
     }

     double solution(double x) 
     {
	  return exp(-x*x);
     }
     
     int test_rk4(int, char**)
     {
	  using namespace cellsim::equation_type;

	  std::cout << "Equation is dy/dx = -2xy" << std::endl;

	  double ynew(1.);
	  double x0 = 0.0;
	  double ynow = 1.0;
	  double da = 0.1;

	  std::cout << "RK4\tAnal. sol\tDeviation" << std::endl;

	  for(unsigned int i(0) ; i < 20 ; ++i)
	  {
	       // std::cout<<"RK4: "<< ynew 
	       // 		<< "\tAnal. sol: "
	       // 		<< solution(x0)
	       // 		<< "\tDev. : "
	       // 		<< fabs(solution(x0) - ynew)
	       // 		<< std::endl;
	       std::cout << std::to_string(ynew)
			 << "\t"
			 << std::to_string(solution(x0))
			 << "\t"
			 << fabs(solution(x0) - ynew)
			 << std::endl;

	       integrators::Standard_RK4<NO_EQN>::do_step(ynow, ynew, x0, da, &f);
	       x0 += da;
	       ynow = ynew;
	  }

	  return 0;
     }
} // namespace test
} // namespace cellsim
