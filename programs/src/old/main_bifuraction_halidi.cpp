#include "main_bifuraction_halidi.hpp"
#include <boost/program_options.hpp>

int biffdiag_halidi(std::string output_pattern,
		    boost::program_options::options_description vm_input)
{
     using cellsim::size_type;
     using namespace cellsim::model_halidi;
     
     constexpr int ncells(1);
     constexpr int n_point(1);
     constexpr double dt(0.01);
     constexpr double dx(1);
     constexpr double t_equ(100);
     constexpr double t_run(1000);
     constexpr double Pip(0.);
     constexpr double Pc(0.);
     constexpr double c_init(0.234333); // < 0.228 mikroM (perhaps down to 100 or 50 nM), 0.225
     constexpr double I_init(0.6); // from 0.8 to 0.9 mikroM, 0.825
     constexpr double s_init(1.28278); // < 1 mikroM (shouldn't go below 0.4mikroM), 0.9
     constexpr double v_init(-45.8178); // from -60 to -44 mV (this is the membrane resting potential), -44
     constexpr double w_init(0.0809169); // < 0.3 (shouldn't go below 0.05), 0.086
     
     using constants_2009::Dc;
     using constants_2009::Dip;

     Cell single_cell(n_point, c_init, I_init, s_init, v_init, w_init, Pip, Pc, Dip, Dc); 
     System mysystem(ncells, single_cell, dt, dx, false);

     std::ofstream cfile("c_output.txt");
     std::ofstream Ifile("I_output.txt");

     std::cout << "Start finding values!" << std::endl;

     using namespace boost::accumulators;
     typedef accumulator_set<double, stats<tag::min, tag::max, tag::mean>> acc_type;

     for(double JPLCago(0.05); JPLCago<=0.2; JPLCago+=0.001){
     	  std::cout << "calculate for JPLCago=" << JPLCago << std::endl;
	  
     	  for(double i(0.) ; i < t_equ ; i += dt){
     	       mysystem.do_step(JPLCago);
     	  }
	  acc_type ca_acc, ip3_acc;

     	  for(double i(0.); i < t_run; i += dt){
     	       mysystem.do_step(JPLCago);
	       ca_acc(mysystem.get_c_at(0));
	       ip3_acc(mysystem.get_I_at(0));
     	  }
	  cfile << JPLCago << " "
		<< mean(ca_acc) << " "
		<< min(ca_acc) << " "
		<< max(ca_acc) << std::endl;
	  Ifile << JPLCago << " "
		<< mean(ip3_acc) << " "
		<< min(ip3_acc) << " "
		<< max(ip3_acc) << std::endl;
     }

     std::cout << "Finished !\n";

     return 0;
}
