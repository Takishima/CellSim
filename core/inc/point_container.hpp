/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POINT_CONTAINER_HPP_INCLUDED
#define POINT_CONTAINER_HPP_INCLUDED

#include "definitions.hpp"
#include "units.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

namespace cellsim {
     
     /*!
      * Container class holding the values of c, I and s for a Cell.
      */
     class Point_Container
     {
     public:
	  typedef std::vector<concentration_ut> container_type;
	  typedef container_type::value_type value_type;
	  typedef container_type::iterator iterator;
	  typedef container_type::const_iterator const_iterator;
	  
	  /*!
	   * Simple constructor
	   *
	   * \param n_point Number of points
	   */
	  Point_Container(size_type n_point);
	  /*!
	   * Simple constructor
	   *
	   * \param n_point Number of points
	   * \param c Value of the cytosolic calcium concentration for all points
	   * \param I Value of the IP3 concentration for all points
	   * \param s Value of the SR calcium concentration for all points
	   */
	  Point_Container(size_type n_point, 
			  value_type c, value_type  I, value_type s);

	  //! Simple accessor on the size of the container
	  size_type size() const;

	  /*!
	   * Access to one particular value
	   *
	   * \tparam type Type of value to get
	   * \param i Index of the point of interest
	   */
	  template <TYPE type>
	  value_type get(size_type i) const;
	  
	  /*!
	   * Swap values of current time with the ones from the next time.
	   *
	   * This function simply swaps the content
	   * of *_array_now_ and *_array_next_
	   */
	  void swap_times();

	  //! Average calcium concentration in cytosol
	  value_type c_avg() const;
	  //! Average IP3 concentration in cytosol
	  value_type I_avg() const;
	  //! Average calcium concentration in the SR
	  value_type s_avg() const;

	  // ====================
	  // First/last element accessors

	  value_type ca_front_now() const;
	  value_type ca_back_now() const;

	  value_type ip3_front_now() const;
	  value_type ip3_back_now() const;

	  value_type s_front_now() const;
	  value_type s_back_now() const;

	  // ====================
	  // Boundaries accessors

	  void set_ca_bc(value_type front, value_type back);
	  void set_ip3_bc(value_type front, value_type back);
	  void set_s_bc(value_type front, value_type back);

	  value_type& ca_bc_front_now();
	  value_type& ca_bc_back_now();
	  value_type& ca_bc_front_next();
	  value_type& ca_bc_back_next();

	  value_type& ip3_bc_front_now();
	  value_type& ip3_bc_back_now();
	  value_type& ip3_bc_front_next();
	  value_type& ip3_bc_back_next();

	  value_type s_bc_front_now() const;
	  value_type s_bc_back_now() const;

	  value_type& s_bc_front_now();
	  value_type& s_bc_back_now();
	  value_type& s_bc_front_next();
	  value_type& s_bc_back_next();


	  // ====================
	  // Iterator accessors
	  const_iterator ca_begin_now() const;
	  const_iterator ca_end_now() const;
	  const_iterator ca_begin_next() const;
	  const_iterator ca_end_next() const;
	  iterator ca_begin_next();
	  iterator ca_end_next();

	  const_iterator ip3_begin_now() const;
	  const_iterator ip3_end_now() const;
	  const_iterator ip3_begin_next() const;
	  const_iterator ip3_end_next() const;
	  iterator ip3_begin_next();
	  iterator ip3_end_next();

	  const_iterator s_begin_now() const;
	  const_iterator s_end_now() const;
	  const_iterator s_begin_next() const;
	  const_iterator s_end_next() const;
	  iterator s_begin_next();
	  iterator s_end_next();

#ifdef TEST
	  iterator begin() {return ca_array_now_.begin();}
	  iterator end() {return ca_array_now_.end();}
#endif // TEST


     private:
	  const size_type last_idx_;

	  container_type ca_array_now_;
	  container_type ca_array_next_;

	  container_type ip3_array_now_;
	  container_type ip3_array_next_;

	  container_type s_array_now_;
	  container_type s_array_next_;
     };
     
     template <>
     inline Point_Container::value_type Point_Container::get<CA>(size_type i) const
     { return ca_array_now_[i]; }
     template <>
     inline Point_Container::value_type Point_Container::get<IP3>(size_type i) const
     { return ip3_array_now_[i]; }
     template <>
     inline Point_Container::value_type Point_Container::get<S>(size_type i) const
     { return s_array_now_[i]; }

} // namespace cellsim

#include "impl/point_container_impl.hpp"

#endif //POINT_CONTAINER_HPP_INCLUDED
