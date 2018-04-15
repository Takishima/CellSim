/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef CONSTANTS_IMPL_HPP_INCLUDED
#define CONSTANTS_IMPL_HPP_INCLUDED

/*!
 * Inspired by BOOST_STATIC_CONSTANT:
 * A convenience macro that allows definition of static
 * constants in headers in an ODR-safe way.
 *
 * \param name Name of the constant
 * \param value Value of the constant
 */
#define STATIC_CONSTANT(name, value)					\
     template<bool b>							\
     struct name##_instance_t						\
     {									\
	  static const std::remove_const<decltype(value)>::type instance; \
     };									\
     									\
     namespace								\
     {									\
	  static const std::remove_const<decltype(value)>::type& name =	\
	       name##_instance_t<true>::instance;			\
     }									\
     									\
     template<bool b>							\
     const std::remove_const<decltype(value)>::type name##_instance_t<b>::instance = value


#endif //CONSTANTS_IMPL_HPP_INCLUDED
