/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef IDENTITY_HPP_INCLUDED
#define IDENTITY_HPP_INCLUDED

namespace cellsim {
namespace integrators {

     /*!
      * Simple template class with an apply method that always returns 
      * its argument.
      */
     template <typename T>
     class Identity
     {
     public:
	  /*!
	   * \param t Value
	   * \return \c t
	   */
	  static T& apply(T& t) {return t;}
     };

     /*!
      * Simple function that represents the identity.
      *
      * \param t Value
      * \return \c t
      */
     template <typename T>
     T& identity(T& t)
	  {
	       return t;
	  }
} // namespace integrators
} // namespace cellsim

#endif //IDENTITY_HPP_INCLUDED
