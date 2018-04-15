#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>

#include "cell_koenigsberger.hpp"
#include "PLC_activator_koenigsberger.hpp"
#include "system_koenigsberger.hpp"

namespace cellsim {
namespace test {

     int test_koenigsberger(int /*argc*/, char** /*argv*/)
     {
	  using namespace model_koenigsberger;
	  
	  size_type Nx(130), Ny(5);
	  time_ut dt(1e-6_s);
	  size_type nstep(20000);
	  size_type nstep_act(5000);
	  
	  // double c(0.1581745);
	  auto c(0.25_uM); // Parthimos et al.
	  auto I(0.1_uM);
	  // double s(0.9);
	  auto s(3_uM); // Parthimos et al.
	  // double v(-64);
	  auto v(-25._mV); // Parthimos et al.
	  // double w(0.0809169);
	  double w(0.5); // Parthimos et al.

	  System system(Nx, Ny, c, I, s, v, w);

	  activator::PLC_Activator activator(0, 10, 0.06_uM / s_);

	  std::ofstream ca_out("ca_out.dat");
	  std::ofstream ip3_out("ip3_out.dat");
	  std::ofstream s_out("s_out.dat");
	  std::ofstream v_out("v_out.dat");

	  std::vector<System::Point> points;
	  points.push_back(std::make_tuple(0, 2));
	  points.push_back(std::make_tuple(4, 2));
	  points.push_back(std::make_tuple(9, 2));
	  points.push_back(std::make_tuple(14, 2));
	  points.push_back(std::make_tuple(19, 2));
	  points.push_back(std::make_tuple(24, 2));
	  points.push_back(std::make_tuple(29, 2));
	  points.push_back(std::make_tuple(100, 2));
	  
	  for(size_type i(0) ; i < nstep ; ++i) {
	       std::cout << "\rEquilibration progress: " << i << "    " << std::flush;
	       system.do_step(dt);
	       system.print_cell_at(points, ca_out, cellsim::CA);
	       system.print_cell_at(points, ip3_out, cellsim::IP3);
	       system.print_cell_at(points, s_out, cellsim::S);
	       system.print_cell_at(points, v_out, cellsim::V);
	  }
	  
	  // for(size_type i(0) ; i < nstep_act ; ++i) {
	  //      std::cout << "\rActivation progress: " << i << "    " << std::flush;
	  //      system.activate(dt, activator);
	  //      system.print_cell_at(points, ca_out, Cell::CA);
	  //      system.print_cell_at(points, ip3_out, Cell::IP3);
	  // }

	  // for(size_type i(0) ; i < nstep ; ++i) {
	  //      std::cout << "\rRelaxation progress: " << i << "    " << std::flush;
	  //      system.do_step(dt);
	  //      system.print_cell_at(points, ca_out, Cell::CA);
	  //      system.print_cell_at(points, ip3_out, Cell::IP3);
	  // }
	  std::cout << "Finished !\n";
	  

	  return 0;
     }
} // namespace test
} // namespace cellsim
