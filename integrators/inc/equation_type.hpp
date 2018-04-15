/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef EQUATION_TYPE_HPP_INCLUDED
#define EQUATION_TYPE_HPP_INCLUDED

#include <type_traits>

#include "definitions.hpp"

namespace cellsim {
namespace equation_type {
     //! No equation
     struct NO_EQN 
     {
	  static constexpr unsigned int value = 100;
     };

     //! No equation
     struct MODEL_PARAM_EQN
     {
	  static constexpr unsigned int value = 101;
     };

     //! All equations
     struct ALL_EQN 
     {
	  static constexpr unsigned int value = 10;
     };

     //! Only the Ca equation
     struct CA_EQN 
     {
	  static constexpr TYPE value = cellsim::CA;
     };
     
     //! Only the IP3 equation
     struct IP3_EQN 
     {
	  static constexpr TYPE value = cellsim::IP3;
     };

     //! Only the s equation
     struct S_EQN 
     {
	  static constexpr TYPE value = cellsim::S;
     };

     //! Only the v equation
     struct V_EQN 
     {
	  static constexpr TYPE value = cellsim::V;
     };

     //! Only the w equation
     struct W_EQN 
     {
	  static constexpr TYPE value = cellsim::W;
     };

     /*!
      * Simple equality metafunction for equation types
      *
      * Base implementation relies on std::is_same
      */
     template <typename U, typename T>
     struct is_same_eq
     {
	  static constexpr bool value = std::is_same<U, T>::value;
     };

     //! Total specialisation for NO_EQN
     template <>
     struct is_same_eq<NO_EQN, NO_EQN>
     { static constexpr bool value = false; };

     //! Partial specialisation for NO_EQN
     template <typename T>
     struct is_same_eq<NO_EQN, T>
     { static constexpr bool value = false; };

     //! Partial specialisation for NO_EQN
     template <typename T>
     struct is_same_eq<T, NO_EQN>
     { static constexpr bool value = false; };

     //! Total specialisation for NO_EQN and ALL_EQN
     template <>
     struct is_same_eq<ALL_EQN, NO_EQN>
     { static constexpr bool value = false; };

     //! Total specialisation for ALL_EQN and NO_EQN
     template <>
     struct is_same_eq<NO_EQN, ALL_EQN>
     { static constexpr bool value = false; };
     
     //! Partial specialisation for ALL_EQN
     template <typename T>
     struct is_same_eq<ALL_EQN, T>
     { static constexpr bool value = true; };

     //! Partial specialisation for ALL_EQN
     template <typename T>
     struct is_same_eq<T, ALL_EQN>
     { static constexpr bool value = true; };

} // namespace integrator

     using equation_type::is_same_eq;

} // namespace cellsim


#endif //EQUATION_TYPE_HPP_INCLUDED
