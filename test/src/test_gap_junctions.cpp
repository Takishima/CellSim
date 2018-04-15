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

int test_gap_junctions(int, char**)
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

     size_type npoints(50);
     double dx(1.);
     double Pip(1.2);
     double Pc(3.25);

     using model_halidi::constants_2009::Dc;
     using model_halidi::constants_2009::Dip;
     
     // _p = prev | _c = current | _n = next
     constexpr double c_p(0.05), I_p(3.0), s_p(1.9);
     constexpr double c_c(0.25), I_c(17.6), s_c(2.4);
     constexpr double c_n(1.15), I_n(20.8), s_n(3.0);

     constexpr double v_p(20);
     constexpr double v_c(23);
     constexpr double v_n(45);
     constexpr double w(0.2);     

     predecessor::Cell cell2009(npoints, c_c, I_c, s_c, v_c, w, Pip, Pc);
     cellsim::model_halidi::Cell cell2012(npoints, c_c, I_c, s_c, v_c, w, Pip, Pc, Dip, Dc);

     double c_front_old, I_front_old, s_front_old, c_back_old, I_back_old, s_back_old;
     double c_front_new, I_front_new, s_front_new, c_back_new, I_back_new, s_back_new;

     std::tie(c_front_old, I_front_old, s_front_old,
	      c_back_old, I_back_old, s_back_old) = cell2009.test_setgap(dx, 
							     cell2009.getboundary_m(), 
							     cell2009.getboundary_p());

     std::tie(c_front_new, I_front_new, s_front_new,
	      c_back_new, I_back_new, s_back_new) = cell2012.test_compute_gap_junctions(dx);

     std::cout << "** Gap junctions computations check WITHOUT neighbours**" << std::endl;

     std::cout << "c_front: "
	       << "\t2009: " << c_front_old
	       << "\tMine: " << c_front_new
	       << "\tDiff: " << std::fabs(c_front_new - c_front_old)
	       << std::endl;
     std::cout << "c_back: "
	       << "\t2009: " << c_back_old
	       << "\tMine: " << c_back_new
	       << "\tDiff: " << std::fabs(c_back_new - c_back_old)
	       << std::endl;

     std::cout << "I_front: "
	       << "\t2009: " << I_front_old
	       << "\tMine: " << I_front_new
	       << "\tDiff: " << std::fabs(I_front_new - I_front_old)
	       << std::endl;
     std::cout << "I_back: "
	       << "\t2009: " << I_back_old
	       << "\tMine: " << I_back_new
	       << "\tDiff: " << std::fabs(I_back_new - I_back_old)
	       << std::endl;

     std::cout << "s_front: "
	       << "\t2009: " << s_front_old
	       << "\tMine: " << s_front_new
	       << "\tDiff: " << std::fabs(s_front_new - s_front_old)
	       << std::endl;
     std::cout << "s_back: "
	       << "\t2009: " << s_back_old
	       << "\tMine: " << s_back_new
	       << "\tDiff: " << std::fabs(s_back_new - s_back_old)
	       << std::endl;


     predecessor::Cell cell2009_p(npoints, c_p, I_p, s_p, v_p, w, Pip, Pc);
     predecessor::Cell cell2009_c(npoints, c_c, I_c, s_c, v_c, w, Pip, Pc);
     predecessor::Cell cell2009_n(npoints, c_n, I_n, s_n, v_n, w, Pip, Pc);

     cellsim::model_halidi::Cell cell2012_p(npoints, c_p, I_p, s_p, v_p, w, Pip, Pc, Dip, Dc);
     cellsim::model_halidi::Cell cell2012_c(npoints, c_c, I_c, s_c, v_c, w, Pip, Pc, Dip, Dc);
     cellsim::model_halidi::Cell cell2012_n(npoints, c_n, I_n, s_n, v_n, w, Pip, Pc, Dip, Dc);

     cell2012_c.set_neighbours(&cell2012_p, &cell2012_n);

     std::tie(c_front_old, I_front_old, s_front_old,
	      c_back_old, I_back_old, s_back_old) = cell2009_c.test_setgap(dx, 
							     cell2009_p.getboundary_m(), 
							     cell2009_n.getboundary_p());

     std::tie(c_front_new, I_front_new, s_front_new,
	      c_back_new, I_back_new, s_back_new) = cell2012_c.test_compute_gap_junctions(dx);

     std::cout << "** Gap junctions computations check with neighbours**" << std::endl;

     std::cout << "c_front: "
	       << "\t2009: " << c_front_old
	       << "\tMine: " << c_front_new
	       << "\tDiff: " << std::fabs(c_front_new - c_front_old)
	       << std::endl;
     std::cout << "c_back: "
	       << "\t2009: " << c_back_old
	       << "\tMine: " << c_back_new
	       << "\tDiff: " << std::fabs(c_back_new - c_back_old)
	       << std::endl;

     std::cout << "I_front: "
	       << "\t2009: " << I_front_old
	       << "\tMine: " << I_front_new
	       << "\tDiff: " << std::fabs(I_front_new - I_front_old)
	       << std::endl;
     std::cout << "I_back: "
	       << "\t2009: " << I_back_old
	       << "\tMine: " << I_back_new
	       << "\tDiff: " << std::fabs(I_back_new - I_back_old)
	       << std::endl;

     std::cout << "s_front: "
	       << "\t2009: " << s_front_old
	       << "\tMine: " << s_front_new
	       << "\tDiff: " << std::fabs(s_front_new - s_front_old)
	       << std::endl;
     std::cout << "s_back: "
	       << "\t2009: " << s_back_old
	       << "\tMine: " << s_back_new
	       << "\tDiff: " << std::fabs(s_back_new - s_back_old)
	       << std::endl;

     
     std::cout.unsetf(std::ios_base::floatfield);
#endif
     return 0;
}

#pragma GCC diagnostic pop

} // namespace test
} // namespace cellsim
