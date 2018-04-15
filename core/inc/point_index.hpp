/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef POINT_INDEX_HPP_INCLUDED
#define POINT_INDEX_HPP_INCLUDED

#include "definitions.hpp"
#include "point_container.hpp"

namespace cellsim {

     /*!
      * \brief Index class to use with Point_Container
      *
      * This class purpose is to provide an index-like functionality to a 
      * Point_Container:
      * \li Access to the elements (via the get method)
      * \li Increment capabilities to navigate across the container
      * \li Get position of index inside the container (via the index method)
      */
     class Point_Index
     {
     public:
	  typedef Point_Container::value_type value_type;

	  /*!
	   * \brief Simple constructor
	   *
	   * \param container Container to link to
	   * \param index_a Starting index
	   * \note Recall that Point_Container stores two extras points
	   *       representing the boundary terms. This is taken into account
	   *       when constructing an Point_Index.
	   */
	  Point_Index(const Point_Container* container, size_type index_a = 0);

	  //! Simple accessor method
	  size_type index() const {return (index_ - 1);}

	  //! Simple accessor method
	  size_type real_index() const {return index_;}

	  /*!
	   * \brief Simple accessor method
	   *
	   * \tparam type Type of value to return
	   * \return Value of type \c type inside the container at the 
	   *         current position
	   */
	  template <TYPE type>
	  value_type get() const {return container_->get<type>(index_);}

	  //! Simple increment operator
	  Point_Index& operator++() {++index_; return *this;}
	  
     private:
	  size_type index_; //!< Position inside the underlying container
	  const Point_Container* container_; //!< Pointer to the actual container
     };
	  
} // namespace cellsim

#endif //POINT_INDEX_HPP_INCLUDED
