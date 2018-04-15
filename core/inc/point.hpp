/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POINT_HPP_INCLUDED
#define POINT_HPP_INCLUDED

#include <cassert>
#include <initializer_list>
#include <iostream>

#include "config.hpp"

namespace cellsim {

     /*!
      * \brief Class representing a discretization point inside a cell.
      *
      * Each cell is divided into a number of Point for the purpose of
      * the simulation.
      */
     class Point
     {
     public:
	  //! Default constructor
	  Point();
	  /*!
	   * Simple constructor
	   *
	   * \param c Calcium concentration in cytosol
	   * \param I IP3 concentration
	   * \param s Calcium concentration in the sarcoplasmic reticulum
	   */
	  Point(double c, double I, double s);
	  
	  //! Simple getter & setter
	  double& c() {return c_;}
	  //! Simple getter
	  double get_c() const {return c_;}
	  //! Simple getter & setter
	  double& I() {return I_;}
	  //! Simple getter
	  double get_I() const {return I_;}
	  //! Simple getter & setter
	  double& s() {return s_;}
	  //! Simple getter
	  double get_s() const {return s_;}

     private:
	  double c_; //!< Calcium concentration in cytosol
	  double I_; //!< IP3 concentration
	  double s_; //!< Calcium concentration in the sarcoplasmic reticulum
     };
}

std::ostream& operator<<(std::ostream& out, cellsim::Point p);

#endif //POINT_HPP_INCLUDED
