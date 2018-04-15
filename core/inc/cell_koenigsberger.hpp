/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Laboratory of Cell Biophysics, 2013"
 * See the LICENSE file for more details.
 */
#ifndef CELL_KOENIGSBERGER_HPP_INCLUDED
#define CELL_KOENIGSBERGER_HPP_INCLUDED

#include "definitions.hpp"
#include "constants_koenigsberger.hpp"
#include "equation_type.hpp"
#include "no_op.hpp"
#include "model_equations_koenigsberger.hpp"
#include "apply_activator.hpp"
#include "units.hpp"

#include <cassert>
#include <iostream>
#include <functional>
#include <tuple>
#include <vector>

namespace cellsim {
     namespace output {
	  class Output_Streams;
     } // namespace output

     namespace model_koenigsberger {
	  
	  class Cell
	  {
	  public:
	       typedef std::vector<const Cell*> neighbours_array;

	       /*!
		* \brief Cell constructor
		*
		* \param c_in Initial calcium concentration in the cytosol
		* \param I_in Initial IP3 concentration
		* \param s_in Initial calcium concentration in the SR
		* \param v_in Initial cell membrane potential
		* \param w_in Initial open state probability of activated potassium channels
		* \param constants Model constants values
		*/
	       Cell(concentration_ut c_in, 
		    concentration_ut I_in,
		    concentration_ut s_in,
		    electric_potential_ut v_in,
		    double w_in,
		    Constants_f_param_t constants);

	       /*!
		* \brief Simple getter/setter of the value of the calcium concentration
		* in the cytosol at the \e current time.
		*/
	       concentration_ut& c() {return c_now_;}
	       /*!
		* \brief Simple getter/setter of the value of the IP3 concentration
		* in the cytosol at the \e current time.
		*/
	       concentration_ut& I() {return I_now_;}
	       /*!
		* \brief Simple getter/setter of the value of the calcium concentration 
		* in the SR in the cytosol at the \e current time.
		*/
	       concentration_ut& s() {return s_now_;}
	       /*!
		* \brief Simple getter/setter of the value of the cell membrane potential
		* in the cytosol at the \e current time.
		*/
	       electric_potential_ut& v() {return v_now_;}
	       /*!
		* \brief Simple getter/setter of the value of the open state probability
		* of activated potassium channels in the cytosol at the \e current time.
		*/
	       double& w() {return w_now_;}

	       //! Constant overload (getter only)
	       concentration_ut c() const {return c_now_;}
	       //! Constant overload (getter only)
	       concentration_ut I() const {return I_now_;}
	       //! Constant overload (getter only)
	       concentration_ut s() const {return s_now_;}
	       //! Constant overload (getter only)
	       electric_potential_ut v() const {return v_now_;}
	       //! Constant overload (getter only)
	       double w() const {return w_now_;}

	       /*!
		* Getter on the neighbours of this Cell.
		*/
	       neighbours_array neighbours() const {return neighbours_;}

	       /*!
		* Simple getter/setter of the value of the calcium concentration
		* in the cytosol at the \e next time.
		*/
	       concentration_ut& c_new() {return c_new_;}
	       /*!
		* Simple getter/setter of the value of the IP3 concentration
		* in the cytosol at the \e next time.
		*/
	       concentration_ut& I_new() {return I_new_;}
	       /*!
		* Simple getter/setter of the value of the calcium concentration 
		* in the SR in the cytosol at the \e next time.
		*/
	       concentration_ut& s_new() {return s_new_;}
	       /*!
		* Simple getter/setter of the value of the cell membrane potential
		* in the cytosol at the \e next time.
		*/
	       electric_potential_ut& v_new() {return v_new_;}
	       /*!
		* Simple getter/setter of the value of the open state probability
		* of activated potassium channels in the cytosol at the \e next time.
		*/
	       double& w_new() {return w_new_;}

	       //! Simple setter
	       void set_neighbours(const neighbours_array array);

	       // //! Simple accessor
	       // Constants_Values& model_constants() {return model_constants_;}
	       // //! Constant overload (getter only)
	       // Constants_Values model_constants() const {return model_constants_;}

	       potential_SR_ut Vcoupling() const;

	       /*!
		* \brief Prints the value of a dynamical variable
		*
		* \param out Output stream
		* \param value Select which dynamical variable to output
		*/
	       void print(std::ostream& out, TYPE value) const;

	       /*!
		* Prints the content of the cell to an output stream
		*
		* \param out Output stream
		*/
	       void print(output::Output_Streams& out) const;

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
		*/
	       void print_fluxes(std::ostream& out,
				 time_ut t) const;


	       /*!
		* Advance one step in time
		*
		* \tparam time_integrator Integrator type to use to integrate
		*                         the simple time derivative  equations.
		* \param t Current time
		* \param dt Time step
		* \param activator Activator functor to use while integrating.\n
		*                 By default, there is no activation (no_op),
		*                 but one can specify any activator here.
		*                 Only the correct equation will be modified
		*                 thanks to the integrator_type type inside
		*                 the activator class.
		*/
	       template <template <typename T> class time_integrator,
			 typename activator_t = activator_impl::no_op<double>>
	       void do_step(time_ut t, 
			    time_ut dt,
			    activator_t& activator);
	       
	       /*!
		* \brief Simple overload in the case of no activators
		*
		*/
	       template <template <typename T> class time_integrator>
	       void do_step(time_ut t, 
			    time_ut dt);

	       /*!
		* Swap quantities at current time (now) with those at future time (new)
		*/
	       void swap_times();

	  private:
	       //! Cytosolic calcium concentration at present time
	       concentration_ut c_now_;
	       //! IP3 concentration at present time
	       concentration_ut I_now_;
	       //! SR calcium concentration at present time
	       concentration_ut s_now_;
	       //! Cell membrane potential at present time
	       electric_potential_ut v_now_;
	       //! Open state probability at present time
	       double w_now_;

	       //! Cytosolic calcium concentration at next time
	       concentration_ut c_new_;
	       //! IP3 concentration at next time
	       concentration_ut I_new_;
	       //! SR calcium concentration at next time
	       concentration_ut s_new_;
	       //! Cell membrane potential at next time
	       electric_potential_ut v_new_;
	       //! Open state probability at next time
	       double w_new_;

	       //! Model constants values
	       Constants_Values model_constants_;
	       
	       //! List of neighbours
	       neighbours_array neighbours_;
	  };
     } // namespace model_keonigsberger
} // namespace cellsim
     
#include "impl/cell_impl.hpp"

template <template <typename T> class time_integrator>
void cellsim::model_koenigsberger::Cell::do_step(time_ut t, 
						 time_ut dt)
{
     static activator_impl::no_op<double> act;
     do_step<time_integrator>(t, dt, act);
}

template <template <typename T> class time_integrator,
	  typename activator_t>
void cellsim::model_koenigsberger::Cell::do_step(time_ut t, 
						 time_ut dt,
						 activator_t& activator)
{
     using namespace equation_type;
     namespace mod_eq = model_koenigsberger::model_equations;
     using namespace std::placeholders;
     using integrators::eval_if;
     typedef eval_if<std::is_same<typename activator_t::equation_type,
				  equation_type::MODEL_PARAM_EQN>::value>
	  model_constants_selector_t;

     const auto& model_constants = model_constants_selector_t::
	  get_model_constants(activator, model_constants_, dt);

     assert(dt.value() > 0.);

     time_integrator<CA_EQN>::do_step(*this, c_new_, 
     				      t, dt, 
				      GET_EQN_PTR(&mod_eq::dc_dt),
     				      //std::bind(&mod_eq::dc_dt, _1, std::cref(model_constants_), _2),
     				      GET_METHOD_PTR_C(&Cell::c),
     				      activator);
     
     time_integrator<IP3_EQN>::do_step(*this, I_new_, 
     				       t, dt, 
     				       GET_EQN_PTR(&mod_eq::dI_dt),
     				       GET_METHOD_PTR_I(&Cell::I),
     				       activator);
     
     time_integrator<S_EQN>::do_step(*this, s_new_, 
     				     t, dt, 
     				     GET_EQN_PTR(&mod_eq::ds_dt),
     				     GET_METHOD_PTR_S(&Cell::s),
     				     activator);
     
     time_integrator<V_EQN>::do_step(*this, v_new_, 
     				     t, dt, 
     				     GET_EQN_PTR(&mod_eq::dv_dt),
     				     GET_METHOD_PTR_V(&Cell::v),
     				     activator);
     
     time_integrator<W_EQN>::do_step(*this, w_new_, 
     				     t, dt,
     				     GET_EQN_PTR(&mod_eq::dw_dt),
     				     GET_METHOD_PTR_W(&Cell::w),
     				     activator);

     // This has to be done when everything has been taken care of :
     //     - integration for this cell
     //     - integration for the rest of the cells.
     //swap_times();
}

#endif //CELL_KOENIGSBERGER_HPP_INCLUDED
