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


void compare_cells2(predecessor::Cell cell2009, cellsim::model_halidi::Cell cell2012, std::string str)
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

int test_model(int, char**)
{
     using namespace model_halidi::model_equations_2009;
     
     const auto c(0.25_uM), I(3.1415_uM), s(2.4_uM);
     Point_Container container(1, c, I, s);
     Point_Index p(&container);
     
     const auto v(23_mV), v_prev(20_mV), v_next(25_mV);
     const auto w(0.2);
     const auto t(0._s);
     const auto JPLCago(flux_ut::from_value(0.06));

     std::cerr << "********************* WARNING *********************\n"
	       << "** Test broken due to constants value incorrect  **\n"
	       << "********************* WARNING *********************\n";

#if 0

     model_halidi::Constants_Values C;

     model_halidi::Cell cell(1, c, I, s, v, w, C);
     
     std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);

     std::cout << "** Model 2009 equations comparison **" << std::endl;

     double old_val(J1(p.get<CA>().value(), p.get<IP3>().value(), p.get<S>().value(), v.value()));
     double new_val(dc_dt(p, v, C, t).value());
     std::cout << "Comparing Ca equation:\n"
	       << "\t2009: " << old_val
	       << "\tMine: " << new_val
	       << "\tDiff: " << std::fabs(old_val - new_val)
	       << std::endl;

     old_val = J2(p.get<CA>().value(), p.get<IP3>().value(), JPLCago.value());
     new_val = dI_dt(p, JPLCago, C, t).value();
     std::cout << "Comparing IP3 equation:\n"
	       << std::scientific
	       << "\t2009: " << old_val
	       << "\tMine: " << new_val
	       << "\tDiff: " << std::fabs(old_val - new_val)
	       << std::endl;
     
     old_val = J3(p.get<CA>().value(), p.get<S>().value());
     new_val = ds_dt(p, C, t).value();
     std::cout << "Comparing s equation:\n"
	       << "\t2009: " << old_val
	       << "\tMine: " << new_val
	       << "\tDiff: " << std::fabs(old_val - new_val)
	       << std::endl;
     
     old_val = J4(p.get<CA>().value(), v.value(), w, v_prev.value(), v_next.value());
     new_val = dv_dt(cell, C, t).value();
     std::cout << "Comparing v equation:\n"
	       << "\t2009: " << old_val
	       << "\tMine: " << new_val
	       << "\tDiff: " << std::fabs(old_val - new_val)
	       << std::endl;

     old_val = J5(p.get<CA>().value(), v.value(), w);
     new_val = dw_dt(cell, C, t).value();
     std::cout << "Comparing w equation:\n"
	       << "\t2009: " << old_val
	       << "\tMine: " << new_val
	       << "\tDiff: " << std::fabs(old_val - new_val)
	       << std::endl << std::endl;
     
     // END of testing model equations
     // ========================================================================
     // BEGIN testing of integration
     
     std::cout << "** Pre-integration comparison **" << std::endl;

     size_type npoints(50);
     auto dt(0.0005_s);
     auto dx(1._um);
     double Pip(1.2);
     double Pc(0.);

     using equation_type::CA_EQN;
     using equation_type::IP3_EQN;
     // using model_halidi::constants_2009::Dc;
     // using model_halidi::constants_2009::Dip;
     
     predecessor::Cell cell2009(npoints, c.value(), I.value(), s.value(), v.value(), w, Pip, Pc);
     cellsim::model_halidi::Cell cell2012(npoints, c, I, s, v, w, C);
     
     // if (cell2009.size() != cell2012.size()) {
     // 	  std::cerr << "Size error: 2009 = " << cell2009.size()
     // 		    << " vs 2012 = " << cell2012.size() << "\n";
     // 	  return 1;
     // }
     // else {
     // 	  std::cout << "Size is " << cell2012.size() << std::endl;
     // }

     compare_cells2(cell2009, cell2012, "");

     std::cout << "** Post-integration comparison **" << std::endl;

     integrators::FTCS<CA_EQN> ca_integrator(dt, dx);
     integrators::FTCS<IP3_EQN> ip3_integrator(dt, dx);

     for (size_type i(0) ; i < 1000 ; ++i) {
	  cell2009.timestep(dt.value(), dx.value(), JPLCago.value());
	  cell2012.do_step<integrators::Forward_Euler>(0._s, dt, JPLCago,
						       ca_integrator,
						       ip3_integrator);
	  cell2012.swap_times();
     }

     compare_cells2(cell2009, cell2012, "");

     std::cout.unsetf(std::ios_base::floatfield);

#endif
     return 0;
}

#pragma GCC diagnostic pop

} // namespace test
} // namespace cellsim
