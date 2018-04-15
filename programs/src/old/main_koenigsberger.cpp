#include "main_koenigsberger.hpp"
#include "units_io.hpp"
#include "model_functions_koenigsberger.hpp"

int main(int /*argc*/, char** /*argv*/)
{
     using cellsim::size_type;
     namespace si = boost::units::si;
     using namespace cellsim;
     using namespace cellsim::model_koenigsberger;
     const auto& s_(si::seconds);
	  
     size_type Nx(130), Ny(5);
     time_ut dt(5e-4_s);
     double t_equ(1);
     double t_act(1);
     double t_run(10);
     size_type nstep_equ(t_equ / dt.value());
     size_type nstep_act(t_act / dt.value());
     size_type nstep_run(t_run / dt.value());
     
     // auto c(0.05_uM);
     auto c(0.1581745_uM);
     //auto c(0.25_uM); // Parthimos et al.
     auto I(0.6_uM);
     // auto s(0.1_uM);
     auto s(3_uM); // Parthimos et al.
     auto v(-65._mV);
     // auto v(-25._mV); Parthimos et al.
     double w(0.0809169);
     // double w(0.5); // Parthimos et al.

     System mysystem(Nx, Ny, dt, c, I, s, v, w);
     activator::PLC_Activator activator(0, 10, 0.35_uM / s_);

     std::ofstream ca_out("ca_out.dat");
     std::ofstream I_out("ip3_out.dat");
     std::ofstream s_out("s_out.dat");
     std::ofstream v_out("v_out.dat");

     size_type idx1(5);
     size_type idx2(99);

     std::cout<< "JPLCago :"  << JPLCago() << std::endl;
     std::cout<< "Nx :"  << Nx << std::endl;
     std::cout<< "Ny :"  << Ny << std::endl;
     std::cout<< "dt :"  << dt << std::endl;
     std::cout<< "nstep_equ :"  << nstep_equ << std::endl;
     std::cout<< "nstep_act :"  << nstep_act << std::endl;
     std::cout<< "nstep_run :"  << nstep_run << std::endl;
     std::cout<< "c :"  << c << std::endl;
     std::cout<< "I :"  << I << std::endl;
     std::cout<< "s :"  << s << std::endl;
     std::cout<< "v :"  << v << std::endl;
     std::cout<< "w :"  << w << std::endl;
     std::cout<< "idx1 :"  << idx1 << std::endl;
     std::cout<< "idx2 :"  << idx2 << std::endl;

     std::ofstream flux_out("fluxes_activ.dat");
     std::ofstream flux_out2("fluxes_normal.dat");
     flux_out << "# t JVOCC JNaCa JNaK JSRuptake JCICR Jextrusion Jleak Jdegrad JIP3 JCl Jback JK Kactiv JPLCago Jact\n"
	      << "# 1   2     3     4     5        6       7        8      9     10   11   12  13   14    15      16"
	      << std::endl;
     flux_out2 << "# t JVOCC JNaCa JNaK JSRuptake JCICR Jextrusion Jleak Jdegrad JIP3 JCl Jback JK Kactiv JPLCago Jact\n"
	       << "# 1   2     3     4     5        6       7        8      9     10   11   12  13   14    15      16"
	       << std::endl;

     ca_out.precision(10);
     I_out.precision(10);
     s_out.precision(10);
     v_out.precision(10);
 
     std::vector<System::Point> points;
     points.push_back(std::make_tuple(5, 2));
     points.push_back(std::make_tuple(10, 2));
     points.push_back(std::make_tuple(20, 2));
     points.push_back(std::make_tuple(100, 2));

     ca_out << "# t idx= ";
     I_out << "# t idx= ";
     s_out << "# t idx= ";
     v_out << "# t idx= ";
     for(auto i: points) {
	  ca_out << std::get<0>(i) << " ";
	  I_out << std::get<0>(i) << " ";
	  s_out << std::get<0>(i) << " ";
	  v_out << std::get<0>(i) << " ";
     }
     ca_out << std::endl;
     I_out << std::endl;
     s_out << std::endl;
     v_out << std::endl;
     
     // points.push_back(std::make_tuple(1, 2));
     // points.push_back(std::make_tuple(2, 2));
     // points.push_back(std::make_tuple(3, 2));
     // points.push_back(std::make_tuple(4, 2));
     // points.push_back(std::make_tuple(5, 2));
     // points.push_back(std::make_tuple(6, 2));
     // points.push_back(std::make_tuple(7, 2));
     // points.push_back(std::make_tuple(8, 2));
     // points.push_back(std::make_tuple(9, 2));
     // points.push_back(std::make_tuple(10, 2));
     // points.push_back(std::make_tuple(11, 2));
     // points.push_back(std::make_tuple(12, 2));

     flux_out.precision(12);
     mysystem.print_cell_at(points, ca_out, cellsim::CA);
     mysystem.print_cell_at(points, I_out, cellsim::IP3);
     mysystem.print_cell_at(points, s_out, cellsim::S);
     mysystem.print_cell_at(points, v_out, cellsim::V);
     mysystem.print_fluxes_at(idx1, 2, flux_out);
     mysystem.print_fluxes_at(idx2, 2, flux_out2);


     for(size_type i(0) ; i < nstep_equ ; ++i) {
     	  std::cout << "\rEquilibration: "
     		    << std::setw(10)  << std::setprecision(4)
     		    << (1.*i/nstep_equ) * 100
     		    << std::flush;
     	  mysystem.do_step();
	  if (i % 10 == 0) {
	       mysystem.print_cell_at(points, ca_out, cellsim::CA);
	       mysystem.print_cell_at(points, I_out, cellsim::IP3);
	       mysystem.print_cell_at(points, s_out, cellsim::S);
	       mysystem.print_cell_at(points, v_out, cellsim::V);
	       mysystem.print_fluxes_at(idx1, 2, flux_out);
	       mysystem.print_fluxes_at(idx2, 2, flux_out2);
	  }
     }
     std::cout << "Finished !\n";

     for(size_type i(0) ; i < nstep_act ; ++i) {
     	  std::cout << "\rActivation: "
     		    << std::setw(10)  << std::setprecision(4)
     		    << (1.*i/nstep_act) * 100
     		    << std::flush;
     	  mysystem.activate(activator);
     	  if (i % 10 == 0) {
     	       mysystem.print_cell_at(points, ca_out, cellsim::CA);
     	       mysystem.print_cell_at(points, I_out, cellsim::IP3);
     	       mysystem.print_cell_at(points, s_out, cellsim::S);
     	       mysystem.print_cell_at(points, v_out, cellsim::V);
     	       mysystem.print_fluxes_at(idx1, 2, activator, flux_out);
     	       mysystem.print_fluxes_at(idx2, 2, activator, flux_out2);
     	  }
     }
     std::cout << "Finished !\n";

     for(size_type i(0) ; i < nstep_run ; ++i) {
     	  std::cout << "\rEquilibration: "
     		    << std::setw(10) << std::setprecision(4)
     		    << (1.*i/nstep_run) * 100
     		    << std::flush;
     	  mysystem.do_step();
     	  if (i % 10 == 0) {
     	       mysystem.print_cell_at(points, ca_out, cellsim::CA);
     	       mysystem.print_cell_at(points, I_out, cellsim::IP3);
     	       mysystem.print_cell_at(points, s_out, cellsim::S);
     	       mysystem.print_cell_at(points, v_out, cellsim::V);
     	       mysystem.print_fluxes_at(idx1, 2, flux_out);
     	       mysystem.print_fluxes_at(idx2, 2, flux_out2);
     	  }
     }
     std::cout << "Finished !\n";


     return 0;
}
