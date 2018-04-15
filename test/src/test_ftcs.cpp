#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>

#include "ftcs.hpp"
#include "point_container.hpp"
#include "point_index.hpp"

namespace cellsim {
namespace test {

     /*
      * u(x, t) = exp(- t) * sin(x)
      *
      * du/dt = D d2u/dx2 + g(u)
      * with w = Dk^2, D = 1 and g(u) = 0
      *
      * u(x,0) = sin(x)
      * u(x,t) = 0 & u(pi/2,t) = 0
      */

     flux_ut g(Point_Index /*i*/, time_ut /*t*/)
     {
	  return flux_ut::from_value(0.);
     }

     double& read(double& x)
     {
	  return x;
     }

     double init(double x)
     {
	  return sin(x);
     }

     double solution(double x, double t)
     {
	  return exp(-t) * sin(x);
     }
     
     int test_ftcs(int, char**)
     {
	  using namespace equation_type;

	  const unsigned int N(200), N1(N+1);
	  auto t(time_ut::from_value(0.));
	  const auto dt(time_ut::from_value(0.00005));
	  const auto dx(length_ut::from_value(M_PI / N1));
	  const auto D(diff_coef_ut::from_value(1.));
	  const auto alpha(D * dt / (dx * dx));
	  const auto factor(1. - 2.*alpha);
	       

	  std::cout << "dx = " << dx << "  dt = " << dt << "  alpha = " << alpha << std::endl;

	  std::vector<double> v, vnew, sol;
	  std::vector<concentration_ut> vnew2;
	  
	  for(unsigned int i(0) ; i <= N1 ; ++i) {
	       v.push_back(init(i*dx.value()));
	  }
	  v.back() = 0.;
	  vnew = v;
	  for (auto el : v) {
	       vnew2.push_back(concentration_ut::from_value(el));
	  }

	  for(unsigned int i(0) ; i <= N1 ; ++i) {
	       sol.push_back(solution(i*dx.value(), dt.value()));
 	  }

	  for (unsigned int i(1) ; i < vnew.size()-1 ; ++i) {
	       vnew[i] = alpha * v.at(i-1) + factor *  v.at(i) + alpha.value() * v.at(i+1);
	  }

	  integrators::FTCS<CA_EQN> integrator(D, dt, dx);

	  Point_Container container(N, 0._uM, 0._uM, 0._uM);

	  if (v.size() != container.size()) {
	       std::cerr << "Size error: v = " << v.size() << " vs container = " << container.size() << "\n";
	       return 1;
	  }
	  
	  auto dest = container.begin();
	  for (auto it(begin(v)) ; it != end(v) ; ++it, ++dest) {
	       *dest = concentration_ut::from_value(*it);
	  }
	  // std::copy(std::begin(v), std::end(v), container.begin());
	  container.ca_bc_front_next() = container.ca_bc_front_now();
	  container.ca_bc_back_next() = container.ca_bc_back_now();

	  integrator.do_step(
	       t, 
	       container,
	       begin<CA, NEXT>(container),
	       &g);

	  vnew2.front() = container.ca_bc_front_next();
	  vnew2.back() = container.ca_bc_back_next();
	  std::copy(begin<CA, NEXT>(container), end<CA, NEXT>(container),
	  	    begin(vnew2)+1);

	  t += dt;

	  std::ofstream out("test_init.dat");
	  for(unsigned int i(0) ; i <= N ; ++i) {
	       out << i*dx.value()
		   << " " << v.at(i) 
		   << " " << vnew.at(i)
		   << " " << vnew2.at(i)
		   << " " << sol.at(i)
		   << std::endl;
	  }

	  using namespace boost::accumulators;
	  typedef accumulator_set<double, stats<tag::min, tag::max>> acc_t;
	  acc_t acc, acc2;

	  std::ofstream diff_out("test_diff.dat");
	  for(unsigned int i(0) ; i <= N1 ; ++i) {
	       diff_out << i << " " 
			<< fabs(sol[i] - vnew[i]) << " "
			<< fabs(sol[i] - vnew2[i].value()) << " "
			<< fabs(vnew[i] - vnew2[i].value()) << std::endl;
	       acc(fabs(sol[i] - vnew[i]));
	       acc2(fabs(sol[i] - vnew2[i].value()));
	  }
	  std::cout << "min = " << min(acc) << "  max = " << max(acc) << std::endl;
	  std::cout << "min2 = " << min(acc2) << "  max2 = " << max(acc2) << std::endl;
	  std::cout << "delta_min = " << fabs(min(acc)-min(acc2))
		    << "  delta_max = " << fabs(max(acc)-max(acc2)) << std::endl;

	  acc = acc_t();

	  return 0;
     }
} // namespace test
} // namespace cellsim
