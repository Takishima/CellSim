/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef CELL_HPP_INCLUDED
#define CELL_HPP_INCLUDED

#include "definitions.hpp"
#include "constants_halidi.hpp"
#include "equation_type.hpp"
#include "apply_activator.hpp"
#include "model_equations_2009.hpp"
#include "model_equations_2012.hpp"
#include "no_op.hpp"
#include "point_container.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <vector>

#ifdef TEST
#include <tuple>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#endif // TEST

GCC_DIAG_OFF(effc++)

namespace cellsim {
     namespace output {
	  class Output_Streams;
	  class Flux_Streams;
     } // namespace output

namespace model_halidi {
     class Cell {
     public:
	  typedef Point_Container container_type;
	  typedef container_type::value_type concentration_t;
	  typedef container_type::iterator iterator;
	  typedef container_type::const_iterator const_iterator;
	  
	  /*!
	   * Cell constructor in the case with no diffusion between cells.
	   *
	   * \param n_point Number of points considered inside the cell
	   *                excluding the boundaries
	   * \param constants Model constants values
	   */
	  Cell(size_type n_point, Constants_f_param_t constants);
	  /*!
	   * Cell constructor in the case with no diffusion between cells.
	   *
	   * \param n_point Number of points considered inside the cell
	   *                excluding the boundaries
	   * \param c_in Initial calcium concentration in the cytosol
	   * \param I_in Initial IP3 concentration
	   * \param s_in Initial calcium concentration in the SR
	   * \param v_in Initial cell membrane potential
	   * \param w_in Initial open state probability of activated potassium 
	   *             channels
	   * \param constants Model constants values
	   */
	  Cell(size_type n_point, 
	       concentration_t c_in,
	       concentration_t I_in,
	       concentration_t s_in,
	       electric_potential_ut v_in, 
	       double w_in, 
	       Constants_f_param_t constants);
	  /*!
	   * \deprecated
	   * \brief Cell constructor.
	   *
	   * \param n_point Number of points considered inside the cell
	   *                excluding the boundaries
	   * \param c_in Initial calcium concentration in the cytosol
	   * \param I_in Initial IP3 concentration
	   * \param s_in Initial calcium concentration in the SR
	   * \param v_in Initial cell membrane potential
	   * \param w_in Initial open state probability of activated potassium channels
	   * \param Pip_in Junctional permeability for IP3
	   * \param Pc_in Junctional permeability for \f$Ca^{2+}\f$
	   * \param Dip_in Diffusion coefficient for IP3
	   * \param Dc_in Diffusion coefficient for \f$Ca^{2+}\f$
	   * \param constants Model constants values
	   */
	  Cell(size_type n_point, 
	       concentration_t c_in,
	       concentration_t I_in,
	       concentration_t s_in,
	       electric_potential_ut v_in, 
	       double w_in, 
	       permeablity_ut Pip_in, 
	       permeablity_ut Pc_in,
	       diff_coef_ut Dip_in,
	       diff_coef_ut Dc_in,
	       Constants_f_param_t constants);

	  //! Average calcium concentration in cytosol inside the cell.
	  concentration_t c_avg() const {return cell_points_.c_avg();}
	  //! Average IP3 concentration in cytosol inside the cell.
	  concentration_t I_avg() const {return cell_points_.I_avg();}
	  //! Average calcium concentration in the SR inside the cell.
	  concentration_t s_avg() const {return cell_points_.s_avg();}

	  //! Simple accessor
	  electric_potential_ut& v() {return v_now_;}
	  electric_potential_ut& v_new() {return v_new_;}
	  //! Simple accessor
	  double& w() {return w_now_;}

	  //! Constant overload (getter only)
	  electric_potential_ut v() const {return v_now_;}
	  //! Constant overload (getter only)
	  double w() const {return w_now_;}

	  // //! Simple accessor
	  // Constants_Values& model_constants() {return model_constants_;}
	  // //! Constant overload (getter only)
	  // Constants_Values model_constants() const {return model_constants_;}

	  /*!
	   * Set the nearest neighbours of this cell
	   *
	   * \param prev The cell just before this one
	   * \param next The cell right after this one
	   */
	  void set_neighbours(const Cell* prev, const Cell* next);

	  /*!
	   * \brief Compute the electrical couplings between cells
	   *
	   * Computes the electrical couplings between this cell and its
	   * neighbours
	   *
	   * \return Value of the electrical couplings
	   */
	  potential_SR_ut Vcoupling() const;

	  /*!
	   * Calculate Ca and IP3 diffusion at gap junctions from previous and 
	   * next cell.\n
	   * Updates the first and last points of the current cell.
	   *
	   * See also report for the exact formula used to take the discontinuity
	   * into account.
	   *
	   * \param dx Spatial step (in micro metres)
	   */
	  void compute_gap_junctions(length_ut dx);

	  /*!
	   * Advance one step in time
	   *
	   * \tparam time_integrator Integrator type to use to integrate
	   *                         the simple time derivative  equations.
	   * \param t Current time
	   * \param dt Time step
	   * \param JPLCago PLC agonist concentration
	   * \param ca_integrator Integrator type to use to integrate
	   *                      the calcium diffusion equation.
	   * \param ip3_integrator Integrator type to use to integrate
	   *                       the IP3 diffusion equation.
	   * \param activator Activator functor to use while integrating.\n
	   *                 By default, there is no activation (no_op),
	   *                 but one can specify any activator here.
	   *                 Only the correct equation will be modified
	   *                 thanks to the integrator_type type inside
	   *                 the activator class.
	   */
	  template <template <typename T> class time_integrator,
		    typename ca_integrator_t,
		    typename ip3_integrator_t,
		    typename activator_type>
	  void do_step(time_ut t, 
		       time_ut dt,
		       flux_ut JPLCago,
		       ca_integrator_t ca_integrator,
		       ip3_integrator_t ip3_integrator,
		       activator_type& activator);

	  //! Simple overload
	  template <template <typename T> class time_integrator,
		    typename ca_integrator_t,
		    typename ip3_integrator_t>
	  void do_step(time_ut t, 
		       time_ut dt,
		       flux_ut JPLCago,
		       ca_integrator_t ca_integrator,
		       ip3_integrator_t ip3_integrator);
	  /*!
	   * Swap values of current time with the ones from the next time.
	   *
	   * This method needs to be called \e before starting the integration
	   * of the next timestep. Otherwise, you will be using the same
	   * values of all the dynamical variables.
	   */
	  void swap_times();

#ifdef TEST
	  std::tuple<double, double, double, double> getca() const;
	  std::tuple<double, double, double, double> getI() const;
	  std::tuple<double, double, double, double> gets() const;

	  std::tuple<double, double, double, double, double, double>
	  test_compute_gap_junctions(length_ut dx);

	  potential_SR_ut Vcoupling(Constants_f_param_t C,
				    electric_potential_ut vm,
				    electric_potential_ut vp) const;


	  double getv() const {return v_now_.value();}
	  double getw() const {return w_now_;}

	  size_type size() const {return cell_points_.size();}
#endif // TEST
	  	  
	  /*!
	   * Prints the content of the cell to an output stream
	   *
	   * \param out Output stream
	   * \param t Current time
	   * \param point_offset Offset to add to position of points
	   */
	  void print_cell(output::Output_Streams& out,
			  time_ut t,
			  size_type& point_offset) const;

	  /*!
	   * Prints the content of one point of the cell to an output stream
	   *
	   * \param out Output stream
	   * \param point_idx Index of point to output
	   */
	  void print_point(output::Output_Streams& out,
			   size_type point_idx) const;

	  /*!
	   * \brief Get column header for fluxes
	   *
	   * \return Formatted string with column headers
	   */
	  static std::string get_fluxes_header();
	  /*!
	   * Prints the content of one point of the cell to an output stream
	   *
	   * \param out Output stream
	   * \param t Current time
	   * \param point_idx Index of point to output
	   */
	  void print_fluxes(std::ostream& out,
			    time_ut t,
			    size_type point_idx) const;

     private:
	  // permeablity_ut Pip_; //!< IP3 permeability
	  // permeablity_ut Pc_; //!< Calcium permeablity
	  // diff_coef_ut Dip_; //!< IP3 diffusion coefficient
	  // diff_coef_ut Dc_; //!< Calcium diffusion coefficient

	  //! Number of points inside the cell (exluding boundary)
	  size_type n_point_;

	  /*!
	   * \brief List of considered points inside a Cell.
	   *
	   * First element is the boundary with previous cell (if exists) and last one
	   * is the boundary with the next cell.
	   */
	  container_type cell_points_;

	  //! Cell membrane potential at present time
	  electric_potential_ut v_now_;
	  //! Open state probability at present time
	  double w_now_;

	  //! Cell membrane potential at next time
	  electric_potential_ut v_new_;
	  //! Open state probability at next time
	  double w_new_;

	  const Cell* prev_; //!< Previous cell in chain (if not nullptr)
	  const Cell* next_; //!< Next cell in chain (if not nullptr)

	  //! Model constants values
	  Constants_Values model_constants_;
     };
} // namespace model_halidi
} // namespace cellsim

GCC_DIAG_ON(effc++)


#include "impl/cell_impl.hpp"

template <template <typename T> class time_integrator,
	  typename ca_integrator_t,
	  typename ip3_integrator_t>
void cellsim::model_halidi::Cell::do_step(time_ut t, 
					  time_ut dt, 
					  flux_ut JPLCago,
					  ca_integrator_t ca_integrator,
					  ip3_integrator_t ip3_integrator)
{
     activator_impl::no_op<double> act;
     do_step<time_integrator>(t, dt, JPLCago, ca_integrator, ip3_integrator, act);
}

template <template <typename T> class time_integrator,
	  typename ca_integrator_t,
	  typename ip3_integrator_t,
	  typename activator_type>
void cellsim::model_halidi::Cell::do_step(time_ut t, 
					  time_ut dt, 
					  flux_ut JPLCago,
					  ca_integrator_t ca_integrator,
					  ip3_integrator_t ip3_integrator,
					  activator_type& activator)
{
     typedef std::function<
	  Point_Index::value_type(const Point_Index)> s_mem_func;
     /*
      * Choose here whether we are using the 2009 or 2012 model
      * Note that you also need to make changes inside
      * compute_gap_junctions()
      */
     using namespace model_halidi;

     assert(dt.value() > 0.);

     using namespace std::placeholders;
     using namespace equation_type;
     using integrators::eval_if;
     typedef eval_if<std::is_same<typename activator_type::equation_type,
				  equation_type::MODEL_PARAM_EQN>::value>
	  model_constants_selector_t;

     const auto& model_constants = model_constants_selector_t::
	  get_model_constants(activator, model_constants_, dt);

#ifdef HALIDI_2009
     namespace mod_eq = model_equations_2009;
#elif defined HALIDI_2012
     namespace mod_eq = model_equations_2012;
#else
# error Unknown or unspecified HALIDI_20XX define macro
#endif
     ca_integrator.setup_system(model_constants.Dc);
     ip3_integrator.setup_system(model_constants.Dip);

     ca_integrator.do_step(
	  t, 
	  cell_points_,
	  begin<CA, NEXT>(cell_points_),
	  std::bind(&mod_eq::dc_dt, _1, v_now_, std::cref(model_constants), _2), 
	  activator);

     ip3_integrator.do_step(
	  t,
	  cell_points_,
	  begin<IP3, NEXT>(cell_points_),
	  std::bind(&mod_eq::dI_dt, _1, JPLCago, std::cref(model_constants), _2), 
	  activator);

     /*
      * In the following, _2 is always a placeholder for the time.
      * (see integrator do_step method)
      */
     assert(cell_points_.size() >= 3);
     auto out_it(begin<S, NEXT>(cell_points_));
     const size_type SIZE(cell_points_.size() - 1);
     for (Point_Index i(&cell_points_) ; i.index() < SIZE; ++i, ++out_it) {
	  time_integrator<S_EQN>::do_step(i, *out_it,
					  t, dt, 
					  std::bind(&mod_eq::ds_dt,
						    _1,
						    std::cref(model_constants),
						    _2),
					  s_mem_func(&Point_Index::get<S>),
					  activator);
     }

     time_integrator<V_EQN>::do_step(*this, v_new_,
				     t, dt,
				     std::bind(&mod_eq::dv_dt, 
					       _1,
					       std::cref(model_constants),
					       _2),
				     GET_METHOD_PTR_V(&Cell::v),
				     activator);

     time_integrator<W_EQN>::do_step(*this, w_new_,
				     t, dt,
				     std::bind(&mod_eq::dw_dt, 
					       _1,
					       std::cref(model_constants),
					       _2),
				     GET_METHOD_PTR_W(&Cell::w),
				     activator);

     v_now_ = v_new_;
     w_now_ = w_new_;
	  
     // This has to be done when everything has been taken care of :
     // 	  - integration
     // 	  - compute_gap_junctions
     // cell_points_.swap_times();
}


#endif //CELL_HPP_INCLUDED
