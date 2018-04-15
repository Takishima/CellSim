#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "point_container.hpp"
#include "point_index.hpp"
#include "cell.hpp"
#include "integrators.hpp"
#include "model_equations_2009.hpp"
#include "read_constants_halidi.hpp"

#include "definition.h"
#include "Cell.h"

namespace cellsim {
namespace test {

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"

     int test_model_equations(int, char**)
     {
	  using namespace model_halidi::model_equations_2009;
     
	  {
	       const auto c(0.25_uM), I(3.1415_uM), s(2.4_uM);
	       Point_Container container(1, c, I, s);
	       Point_Index p(&container);
     
	       const auto v(23_mV), v_prev(20_mV), v_next(25_mV);
	       const auto w(0.2);
	       const auto t(0._s);
	       const auto JPLCago(flux_ut::from_value(0.06));

	       model_halidi::Constants_Values C;
	       if (!input::read_constants("halidi_constants.cfg", C)) {
		    std::cerr << "Error while reading constants file\n";
		    return 1;
	       }

	       model_halidi::Cell cell(1, c, I, s, v, w, C);

	       std::cout.precision(14);
	       std::cout << "ca  = " << dc_dt(p, v, C, t).value() << std::endl;
	       std::cout << "ip3 = " << dI_dt(p, JPLCago, C, t).value() << std::endl;
	       std::cout << "s   = " << ds_dt(p, C, t).value() << std::endl;
	       std::cout << "v   = " << dv_dt(cell, C, t).value() << std::endl;
	       std::cout << "v'  = " << dv_dt_test(cell, v_prev, v_next, C, t).value() << std::endl;
	       std::cout << "w   = " << dw_dt(cell, C, t).value() << std::endl;
	       std::cout << std::endl;
	  }

	  {
	       const auto c(0.34_uM), I(2.2384_uM), s(34.4_uM);
	       Point_Container container(1, c, I, s);
	       Point_Index p(&container);
     
	       const auto v(-64._mV), v_prev(-234._mV), v_next(34.5_mV);
	       const auto w(0.3123);
	       const auto t(0._s);
	       const auto JPLCago(flux_ut::from_value(10.23));

	       model_halidi::Constants_Values C;
	       if (!input::read_constants("halidi_constants.cfg", C)) {
		    std::cerr << "Error while reading constants file\n";
		    return 1;
	       }

	       model_halidi::Cell cell(1, c, I, s, v, w, C);

	       std::cout.precision(14);
	       std::cout << "ca  = " << dc_dt(p, v, C, t).value() << std::endl;
	       std::cout << "ip3 = " << dI_dt(p, JPLCago, C, t).value() << std::endl;
	       std::cout << "s   = " << ds_dt(p, C, t).value() << std::endl;
	       std::cout << "v   = " << dv_dt(cell, C, t).value() << std::endl;
	       std::cout << "v'  = " << dv_dt_test(cell, v_prev, v_next, C, t).value() << std::endl;
	       std::cout << "w   = " << dw_dt(cell, C, t).value() << std::endl;
	       std::cout << std::endl;
	  }
	  return 0;
     }
} // namespace test
} // namespace cellsim
