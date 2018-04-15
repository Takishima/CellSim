#if 0
#include <iostream>
#include <fstream>
#include <functional>
#include <tuple>

#include "cell.hpp"
#include "integrators.hpp"
#include "PLC_activator.hpp"
#include "system.hpp"
#include "no_op.hpp"
#include "cell.hpp"
#include "point.hpp"

#include "model_equations_2009.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"

int main()
{     
     using namespace cellsim;
     typedef std::size_t size_type;
     using namespace cellsim::model_halidi;

     size_type ncells(15);
     size_type npoints(50);
     //size_type npoints(500);
     
     // Activation param
     size_type cell_idx(7);
     size_type act_region(10);
     
     auto dt(time_ut::from_value(1.01e-3));		//--time step
     // double dt(0.0005);		//--time step
     auto dx(length_ut::from_value(1.));			//--space step
     auto JPLCago(flux_ut::from_value(0.06));	//--JPLCagonist flux
     double t_equ(0.01);		//--duration for finding the equilibrium position
     size_type n_equ(t_equ/dt.value());
     double t_act(8);		//--duration of the activation of the c.w.
     size_type n_act(t_act/dt.value());
     double t_evo(160);		//--duration of the evolution after the activation
     size_type n_evo(t_evo/dt.value());
     auto Pip(permeablity_ut::from_value(1.1));			//--IP_3 permeability
     auto Pc(permeablity_ut::from_value(10));			//--calcium permeablity
     //double J_act(0.12);		//--activation flux
     auto J_act(flux_ut::from_value(0.06));		//--activation flux
     size_type o_n(200);			//--number that defines the output times

     
     auto c_init(concentration_ut::from_value(0.234333)); // < 0.228 mikroM (perhaps down to 100 or 50 nM),   0.225
     auto I_init(concentration_ut::from_value(0.6));      // from 0.8 to 0.9 mikroM,    0.825
     auto s_init(concentration_ut::from_value(1.28278));  // < 1 mikroM (shouldn't go below 0.4mikroM)   0.9
     auto v_init(electric_potential_ut::from_value(-45.8178)); // from -60 to -44 mV (this is the membrane resting potential),	-44
     auto w_init(0.0809169);// < 0.3 (shouldn't go below 0.05)    0.086					

     using constants_2009::Dip;
     using constants_2009::Dc;

     Cell single_cell(npoints, c_init, I_init, s_init, v_init, w_init, Pip, Pc, Dip, Dc);

     System system(ncells, single_cell, dt, dx, false);


     std::cout << "dt = " << dt << std::endl
	       << "dx = " << dx << std::endl
	       << "JPLCago = " << JPLCago << std::endl
	       << "t_equ = " << t_equ << std::endl
	       << "n_equ = " << n_equ << std::endl
	       << "t_act = " << t_act << std::endl
	       << "t_evo = " << t_evo << std::endl
	       << "Pip = " << Pip << std::endl
	       << "Pc = " << Pc << std::endl
	       << "Dip = " << Dip << std::endl
	       << "Dc = " << Dc << std::endl
	       << "J_act = " << J_act << std::endl
	       << "c_init = " << c_init << std::endl
	       << "I_init = " << I_init << std::endl
	       << "s_init = " << s_init << std::endl
	       << "v_init = " << v_init << std::endl
	       << "w_init = " << w_init << std::endl
	       << "center = " << npoints/2 << std::endl
	  ;

//-----------------open output files-----------------
     std::ofstream cfile("P" + std::to_string(Pip.value()) + "-0_c_output.txt");
     std::ofstream Ifile("P" + std::to_string(Pip.value()) + "-0_I_output.txt");
     std::ofstream vfile("P" + std::to_string(Pip.value()) + "-0_v_output.txt");
     cfile.precision(12);
     Ifile.precision(12);
     vfile.precision(12);

     system.print(cfile, Ifile, vfile);

     for(size_type i(0) ; i < n_equ ; ++i) {
     	  system.do_step(JPLCago);

     	  if (i % o_n == 0){
     	       system.print(cfile, Ifile, vfile);
     	       std::cerr << "\r" << "Looking for equilibrium........" 
     			 << ((dt.value()*i)/t_equ*100) << "%    " << std::flush;
     	  }
     }
     std::cerr << "\r" << "Looking for equilibrium........at equilibrium!"  << std::endl;

     activator::PLC_Activator activator(cell_idx, npoints/2, act_region, J_act);
     for(size_type i(0) ; i < n_act ; ++i) {
	  system.do_step(JPLCago, activator);

	  if (i % o_n == 0){
	       system.print(cfile, Ifile, vfile);
	       std::cerr << "\r" << "Activation in progress........." <<
		    ((dt.value()*i)/t_act*100) << "%    " << std::flush;
	  }
     }
     std::cerr << "\r" << "Activation in progress.........activiation finished!" << std::endl;

     std::cerr << "Evolution in progress........" << std::flush;
     for(size_type i(0) ; i < n_evo ; ++i){
	  system.do_step(JPLCago);

	  if (i % o_n == 0){
	       system.print(cfile, Ifile, vfile);
	       std::cerr << "\r" << "Evolution in progress.........." << 
		    ((dt.value()*i)/t_evo*100) << "%    " << std::flush;
	  }
     }
     std::cerr << "\r" << "Evolution in progress..........program finished!" << std::endl;

     
     return 0;
}
#pragma GCC diagnostic pop

#endif 
