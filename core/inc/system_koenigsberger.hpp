/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef SYSTEM_KOENIGSBERGER_HPP_INCLUDED
#define SYSTEM_KOENIGSBERGER_HPP_INCLUDED

#include "config.hpp"
#include "cell_koenigsberger.hpp"
#include "integrators.hpp"
#include "units.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#ifdef MULTI_THREAD
#  include <thread>
#endif // MULTI_THREAD

#ifdef WITH_RANGE_CHECK
#  include <stdexcept>
#endif // WITH_RANGE_CHECK

namespace cellsim {
     namespace output {
	  class Output_Streams;
	  class Flux_Streams;
     } // namespace output

     namespace model_koenigsberger {	  
	  /*!
	   * \brief Class representing the system of cells we are simulating.
	   *
	   * Internally, the 2D arrangment of cells is stored in a 1D array:
	   * \verbatim
	    0    1    2    3   ... Nx-1
	    Nx  Nx+1 Nx+2 Nx+3 ... 2Nx-1
	    ...
	    \endverbatim
	   */
	  class System
	  {
	  public:
	       typedef std::vector<Cell> cell_array;
	       typedef cell_array::iterator cell_iterator;
	       // typedef std::unique_ptr<std::ofstream> ofstream_ptr;
	       // typedef std::pair<ofstream_ptr, Cell::VALUE> IO_pair;
	       typedef std::tuple<size_type, size_type> Point;

	       /*!
		* \brief Simple constructor
		*
		* \param Nx Number of cells in the x-direction
		* \param Ny Number of cells in the y-direction
		* \param constants Model constants values
		*/
	       System(size_type Nx, size_type Ny, Constants_f_param_t constants);
	       /*!
		* \brief Simple constructor
		*
		* \param Nx Number of cells in the x-direction
		* \param Ny Number of cells in the y-direction
		* \param c Initial calcium concentration in the cytosol
		* \param I Initial IP3 concentration
		* \param s Initial calcium concentration in the SR
		* \param v Initial cell membrane potential
		* \param w Initial open state probability of activated potassium channels
		* \param constants Model constants values
		*/
	       System(size_type Nx, 
		      size_type Ny,
		      concentration_ut c, 
		      concentration_ut I,
		      concentration_ut s, 
		      electric_potential_ut v, 
		      double w, 
		      Constants_f_param_t constants);

	       /*!
		* \brief Advance one step in time
		*/
	       void do_step(time_ut dt);
	       
	       /*!
		* \brief Advance one step in time with activation
		*
		* Advance one step in time while activating some of the cells.
		* (will be removed in a future revision in favour of do_step)
		*/
	       template <typename time_t, typename activator_t>
	       void activate(time_t dt, activator_t& activator);

	       /*!
		* \brief Print some information about a cell
		*
		* \param x X coordinate of the cell to print
		* \param y Y coordinate of the cell to print
		* \param out Output stream
		* \param value Which dynamical variable to print out
		*/
	       void print_cell_at(size_type x, size_type y,
				  std::ofstream& out, TYPE value) const;
	       
	       //! Simple overload
	       void print_cell_at(size_type x, size_type y,
				  output::Output_Streams& out) const;

	       /*!
		* \brief Simple overload
		*
		* Allows to print a bunch of cells at once.
		* \param array Array of Points containing the coordinates of
		*              the cells to print
		* \param out Output stream
		* \param value Which dynamical variable to print out
		*/
	       void print_cell_at(std::vector<Point> array,
				  std::ofstream& out, TYPE value) const;
	       
	       //! Simple overload
	       void print_cell_at(std::vector<Point> array,
				  output::Output_Streams& out) const;


	       /*!
		* \brief Print the values of the different fluxes
		*
		* Prints the value of all the fluxes for a particular cell
		*
		* \param p Indices of the cell to output data for
		* \param out Output stream
		*/
	       void print_fluxes_at(Point p,
				    std::ofstream& out) const;
	       
	       /*!
		* \brief Simple overload
		*/
	       void print_fluxes_at(std::vector<Point> array,
				    output::Flux_Streams& out) const;
	       /*!
		* \brief Accessor to calcium concentration
		*
		* \param x X-index of a cell
		* \param y Y-index of a cell
		* \return Calcium concentration inside cell
		*/
	       concentration_ut get_c_at(size_type x, size_type y) const;
	       /*!
		* \brief Simple overload for IP3
		*
		* \param x X-index of a cell
		* \param y Y-index of a cell
		* \return IP3 concentration inside cell
		*/
	       concentration_ut get_I_at(size_type x, size_type y) const;

	       // void print(std::ofstream& out, TYPE value) const;
	       // void print(std::vector<IO_pair>& array) const;

	  private:
	       //! Simple helper function
	       void print_fluxes_helper_(size_type x, size_type y,
					 std::ofstream& out) const;
	       
	       /*!
		* \brief Method to call after any time stepping
		*/
	       void end_step_(time_ut dt);

	       /*!
		* \brief Cell accessor
		*
		* \param x X coordinate of the cell to print
		* \param y Y coordinate of the cell to print
		*/
	       Cell& cell_at_ (size_type x, size_type y);
	       /*!
		* \brief Const overload
		*
		* \param x X coordinate of the cell to print
		* \param y Y coordinate of the cell to print
		*/
	       Cell cell_at_ (size_type x, size_type y) const;

	       /*!
		* \brief Cell address accessor
		*
		* \param x X coordinate of the cell to print
		* \param y Y coordinate of the cell to print
		*/
	       const Cell* address_at_(size_type x, size_type y) const;

	       /*!
		* \brief Set neighbours of all cells
		*/
	       void set_neighbours_();

#ifdef WITH_RANGE_CHECK
	       //! Checks whether \c x and \c y are out of range or not
	       void range_check_(size_type x, size_type y) const;
#endif // WITH_RANGE_CHECK

	       const size_type Nx_; //!< Number of cells in the X-direction
	       const size_type Ny_; //!< Number of cells in the Y-direction

	       time_ut t_; //!< Current time
	       
	       cell_array cells_; //!< Container of cells
	  };
     } // namespace model_keonigsberger
} // namespace cellsim

#include "impl/system_koenigsberger_impl.hpp"

#endif //SYSTEM_KOENIGSBERGER_HPP_INCLUDED
