/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef APPLY_ACTIVATOR_HPP_INCLUDED
#define APPLY_ACTIVATOR_HPP_INCLUDED

namespace cellsim {
namespace integrators {
    /*!
      * Template metafunction to use or not use the activator
      *
      * This class basically let us choose whether we call the 
      * activator function method or not based on the value of
      * \c use_activator at \b compile \b time.
      */
     template <bool use_activator>
     struct eval_if
     {
	  template <typename function_t, 
		    typename value_t,
		    typename variable_t,
		    typename activator_t>
	  static auto apply_activator(const function_t f,
				      const value_t v,
				      const variable_t a,
				      activator_t& /*act_func*/)
	       -> decltype(f(v,a))
	       {
		    return f(v, a);
	       }

	  //! Special case for NO_EQN & MODEL_PARAM_EQN type activators
	  template <typename activator_t,
		    typename value_t,
		    typename time_t>
	  static void apply_activator(activator_t& /*act_func*/,
				      value_t& /*v*/,
				      time_t /*dt*/)
	  {}

	  template <typename activator_t,
		    typename constants_t,
		    typename time_t>
	  static constants_t get_model_constants(activator_t& /*a*/,
						 constants_t& c,
						 time_t /*dt*/)
	       {
		    return c;
	       }
     };

     /*!
      * Template specialisation in the case we \e do want to use the 
      * activator.
      */
     template <>
     struct eval_if<true>
     {
	  template <typename function_t, 
		    typename value_t,
		    typename variable_t,
		    typename activator_t>
	  static auto apply_activator(const function_t f,
				      const value_t v,
				      const variable_t a,
				      activator_t& act_func)
	       -> decltype(f(v,a))
	       {
		    return act_func(f, v, a);
	       }

	  //! Special case for NO_EQN & MODEL_PARAM_EQN type activators
	  template <typename activator_t,
		    typename value_t,
		    typename time_t>
	  static void apply_activator(activator_t& act_func,
				      value_t& v,
				      time_t dt)
	       {
		    act_func(v, dt);
	       }

	  template <typename activator_t,
		    typename constants_t,
		    typename time_t>
	  static constants_t get_model_constants(activator_t& a,
						 constants_t& /*c*/,
						 time_t dt)
	       {
		    return a.get_constants(dt);
	       }
     };
} // namespace integrators
} // namespace cellsim

#endif //APPLY_ACTIVATOR_HPP_INCLUDED
