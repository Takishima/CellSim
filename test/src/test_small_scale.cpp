#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "point_container.hpp"
#include "point_index.hpp"
#include "cell.hpp"
#include "integrators.hpp"
#include "model_equations_2009.hpp"

#include "definition.h"
#include "Cell.h"

namespace cellsim {
namespace test {

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"


void compare_cells(predecessor::Cell cell2009, cellsim::model_halidi::Cell cell2012, std::string str)
{
     double std_old(0.), std_new(0.), min_old(0.), min_new(0.), max_old(0.), max_new(0.);
     double old_val(0.), new_val(0.);

     std::tie(old_val, std_old, min_old, max_old) = cell2009.getca();
     std::tie(new_val, std_new, min_new, max_new) = cell2012.getca();
     std::cout << "Ca " << str << ":\n"
	       << "\t2009: " << old_val
	       << "  2012: " << new_val
	       // << "  Diff: " << std::fabs(old_val - new_val)
	       << "  min_diff: " << std::fabs(min_old - min_new)
	       << "  max_diff: " << std::fabs(max_old - max_new)
	       << std::endl;

     std::tie(old_val, std_old, min_old, max_old) = cell2009.getI();
     std::tie(new_val, std_new, min_new, max_new) = cell2012.getI();
     std::cout << "IP3 " << str << ":\n"
	       << "\t2009: " << old_val
	       << "  2012: " << new_val
	       // << "  Diff: " << std::fabs(old_val - new_val)
	       << "  min_diff: " << std::fabs(min_old - min_new)
	       << "  max_diff: " << std::fabs(max_old - max_new)
	       << std::endl;

     std::tie(old_val, std_old, min_old, max_old) = cell2009.gets();
     std::tie(new_val, std_new, min_new, max_new) = cell2012.gets();
     std::cout << "s " << str << ":\n"
	       << "\t2009: " << old_val
	       << "  2012: " << new_val
	       // << "  Diff: " << std::fabs(old_val - new_val)
	       << "  min_diff: " << std::fabs(min_old - min_new)
	       << "  max_diff: " << std::fabs(max_old - max_new)
	       << std::endl;

     old_val = cell2009.getv();
     new_val = cell2012.getv();
     std::cout << "v " << str << ":\n"
	       << "\t2009: " << old_val
	       << "  2012: " << new_val
	       << "  Diff: " << std::fabs(old_val - new_val)
	       << std::endl;

     old_val = cell2009.getw();
     new_val = cell2012.getw();
     std::cout << "w " << str << ":\n"
	       << "\t2009: " << old_val
	       << "  2012: " << new_val
	       << "  Diff: " << std::fabs(old_val - new_val)
	       << std::endl;

     std::cout << std::endl;
}

int test_small_scale(int, char**)
{
     /*
      *
      *
      * TEST BROKEN
      *
      *
      */
     
     
     std::cerr << "********************* WARNING *********************\n"
	       << "** Test broken due to constants value incorrect  **\n"
	       << "********************* WARNING *********************\n";

#if 0
     using namespace model_halidi::model_equations_2009;

     std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
     std::cout.precision(10);

     const size_type nstep(1000);
     size_type npoints(50);
     double dx(1.);
     double dt(0.0005);
     double JPLCago(0.06);
     double Pip(1.2);
     double Pc(3.25);

     // using model_halidi::constants_2009::Dc;
     // using model_halidi::constants_2009::Dip;
     
     // _p = prev | _c = current | _n = next
     constexpr double c_p(0.05), I_p(3.0), s_p(1.9);
     constexpr double c_c(0.25), I_c(17.6), s_c(2.4);
     constexpr double c_n(1.15), I_n(20.8), s_n(3.0);

     constexpr double v_p(20);
     constexpr double v_c(23);
     constexpr double v_n(45);
     constexpr double w(0.2); 

     using equation_type::CA_EQN;
     using equation_type::IP3_EQN;
     integrators::FTCS<CA_EQN> ca_integrator(dt, dx);
     integrators::FTCS<IP3_EQN> ip3_integrator(dt, dx);    

     predecessor::Cell cell2009_p(npoints, c_p, I_p, s_p, v_p, w, Pip, Pc);
     predecessor::Cell cell2009_c(npoints, c_c, I_c, s_c, v_c, w, Pip, Pc);
     predecessor::Cell cell2009_n(npoints, c_n, I_n, s_n, v_n, w, Pip, Pc);

     cellsim::model_halidi::Cell cell2012_p(npoints, c_p, I_p, s_p, v_p, w, Pip, Pc, Dip, Dc);
     cellsim::model_halidi::Cell cell2012_c(npoints, c_c, I_c, s_c, v_c, w, Pip, Pc, Dip, Dc);
     cellsim::model_halidi::Cell cell2012_n(npoints, c_n, I_n, s_n, v_n, w, Pip, Pc, Dip, Dc);

     cell2012_p.set_neighbours(nullptr, &cell2012_c);
     cell2012_c.set_neighbours(&cell2012_p, &cell2012_n);
     cell2012_n.set_neighbours(&cell2012_c, nullptr);

     std::cout << "** Pre-integration comparison **" << std::endl;

     compare_cells(cell2009_p, cell2012_p, "previous");
     compare_cells(cell2009_c, cell2012_c, "current");
     compare_cells(cell2009_n, cell2012_n, "next");

     // for(int j=0; j<=n_cell; ++j){
     // 	  l[j].timestep(dt, dx, JPLCago);
     // }

     // l[0].setgap(dx, l[0].getboundary_m(), l[1].getboundary_m());
     // for(int j=1; j<=n_cell-1; ++j){
     // 	  l[j].setgap(dx, l[j-1].getboundary_p(), l[j+1].getboundary_m());
     // };
     // l[n_cell].setgap(dx, l[n_cell-1].getboundary_p(), l[n_cell].getboundary_p());
     
     std::cout << "****    Integrating for "
	       << nstep << " steps    ****" << std::endl << std::endl;

     for (unsigned int i(0) ; i < nstep ; ++i) {
	  // Stepping 2009 cells
	  cell2009_p.timestep(dt, dx, JPLCago);
	  cell2009_c.timestep(dt, dx, JPLCago);
	  cell2009_n.timestep(dt, dx, JPLCago);

	  // Stepping 2012 cells
     	  cell2012_p.do_step<integrators::Forward_Euler>(0., dt, JPLCago,
							 ca_integrator,
							 ip3_integrator);
	  cell2012_c.do_step<integrators::Forward_Euler>(0., dt, JPLCago,
	  						 ca_integrator,
	  						 ip3_integrator);
	  cell2012_n.do_step<integrators::Forward_Euler>(0., dt, JPLCago,
	  						 ca_integrator,
	  						 ip3_integrator);

	  cell2009_p.setgap(dx, cell2009_p.getboundary_m(), cell2009_c.getboundary_m());
	  cell2009_c.setgap(dx, cell2009_p.getboundary_p(), cell2009_n.getboundary_m());
	  cell2009_n.setgap(dx, cell2009_c.getboundary_p(), cell2009_n.getboundary_p());

	  cell2012_p.swap_times();
	  cell2012_c.swap_times();
	  cell2012_n.swap_times();

	  cell2012_p.compute_gap_junctions(dx);
	  cell2012_c.compute_gap_junctions(dx);
	  cell2012_n.compute_gap_junctions(dx);

     }

     std::cout << "** Post-integration comparison **" << std::endl;

     compare_cells(cell2009_p, cell2012_p, "previous");
     compare_cells(cell2009_c, cell2012_c, "current");
     compare_cells(cell2009_n, cell2012_n, "next");
     

     std::cout.unsetf(std::ios_base::floatfield);
#endif
     return 0;
}

#pragma GCC diagnostic pop

} // namespace test
} // namespace cellsim
