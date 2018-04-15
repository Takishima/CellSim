/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef SYSTEM_HPP_INCLUDED
#define SYSTEM_HPP_INCLUDED

#include "config.hpp"
#include "cell.hpp"
#include "integrators.hpp"
#include "equation_type.hpp"
#include "units.hpp"

#include <cassert>
#include <functional>
#include <stdexcept>
#include <ostream>
#include <vector>

#ifdef WITH_RANGE_CHECK
#  include <stdexcept>
#endif // WITH_RANGE_CHECK

namespace cellsim {
     namespace output {
	  class Output_Streams;
	  class Flux_Streams;
     } // namespace output

namespace model_halidi {
     class System
     {
     public:
	  typedef size_type Point;
	  template <typename T>
	  using integrator_t = integrators::Forward_Euler<T>;

	  /*!
	   * Constructor
	   *
	   * \param Ncells The total number of cells in the system
	   * \param c A cell example (will be copied \c Ncells times)
	   * \param dx_in Integration spatial step
	   * \param periodic_bc Use periodic boundary conditions ?
	   */
	  System(size_type Ncells,
		 Cell c,
		 length_ut dx_in,
		 bool periodic_bc = true);

	  //! Simple accessor
	  length_ut dx() const {return dx_;}

	  /*!
	   * \brief Simple setter
	   *
	   * \param dx_in New spatial step
	   */
	  void set_dx(length_ut dx_in) 
	       { assert(dx_in.value() > 0); dx_ = dx_in; }
	  
	  /*!
	   * Integrate all differential equations one time step.
	   */
	  void do_step(time_ut dt,
		       flux_ut JPLCago);

	  template <typename activator_t>
	  void activate(time_ut dt,
			flux_ut JPLCago,
			activator_t& activator);
	  
	  /*!
	   * Print the content of the system to output streams
	   */
	  void print_system(output::Output_Streams& out) const;
	  void print_cell_at(size_type cell_idx,
			     size_type point_idx,
			     output::Output_Streams& out) const;

	  void print_cell_at(std::vector<Point> array,
			     size_type point_idx,
			     output::Output_Streams& out) const;
	  
	  void print_fluxes_at(size_type cell_idx,
			       size_type point_idx,
			       std::ostream& out) const;

	  void print_fluxes_at(std::vector<Point> array,
			       size_type point_idx,
			       output::Flux_Streams& out) const;

	  concentration_ut get_c_at(size_type x) const;
	  concentration_ut get_I_at(size_type x) const;

	  static std::string get_system_header();
	  
     private:
	  /*!
	   * \brief Cell accessor
	   *
	   * \param i Index of cell
	   */
	  Cell& cell_at_ (size_type i);
	  /*!
	   * \brief Const overload
	   *
	   * \param i Index of cell
	   */
	  Cell cell_at_ (size_type i) const;

	  /*!
	   * Set neighbours for each cell of the system
	   *
	   * \param periodic_bc Use periodic boundary conditions ?
	   */
	  void set_neighbours_(bool periodic_bc);

#ifdef WITH_RANGE_CHECK
	       //! Checks whether \c i is out of range or not
	       void range_check_(size_type i) const;
#endif // WITH_RANGE_CHECK

	  //! Method to call after each time stepping
	  void post_update_step_(time_ut dt);
	  	  
	  time_ut t_; //!< Current time
	  length_ut dx_; //!< Spatial step
	  std::vector<Cell> cells_; //!< Array of cells
     };

     template <typename activator_t>
     void System::activate(time_ut dt,
			   flux_ut JPLCago, 
			   activator_t& activator)
     {
	  using namespace cellsim::equation_type;
	  using integrators::eval_if;
	  using equation_type::is_same_eq;

	  integrators::FTCS<CA_EQN> ca_integrator(dt, dx_);
	  integrators::FTCS<IP3_EQN> ip3_integrator(dt, dx_);

	  /*
	   * First take care of equations for Ca, IP3, v and w
	   */
	  const size_type SIZE(cells_.size());
	  for(size_type i(0) ; i < SIZE ; ++i) {
	       if (activator.do_activation(i)) {
		    cells_[i].do_step<integrator_t>(t_, dt, 
						    JPLCago,
						    ca_integrator,
						    ip3_integrator,
						    activator);
		    /*
		     * If activator is designed to act directly on a cell,
		     * apply it now that integration is done
		     */
		    eval_if<std::is_same<typename activator_t::equation_type,
					 equation_type::NO_EQN>::value>::
			 apply_activator(activator, cells_[i], dt);
	       }
	       else {
		    cells_[i].do_step<integrator_t>(t_, dt, 
						    JPLCago,
						    ca_integrator,
						    ip3_integrator);
	       }
	       cells_[i].compute_gap_junctions(dx_);
	  }
	  // for (auto& c : cells_) {
	  //      // advance to next step
	  //      c.do_step<integrator_t>(t_, dt, JPLCago,
	  // 					     ca_integrator,
	  // 					     ip3_integrator);
	  //      c.swap_times();
	  // }

	  post_update_step_(dt);
     }

} // namespace halidi
} // namespace cellsim

// #include "impl/system_impl.hpp"

#endif //SYSTEM_HPP_INCLUDED
