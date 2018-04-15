/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef UNITS_IMPL_HPP_INCLUDED
#define UNITS_IMPL_HPP_INCLUDED

#ifdef NO_UNITS

#define MAKE_LITERAL_OPERATOR(suffix)					\
     inline units_impl::quantity_t operator"" suffix(const long double a) \
     {									\
	  return units_impl::quantity_t::from_value(a);			\
     }									\
     inline units_impl::quantity_t operator"" suffix(const unsigned long long a) \
     {									\
	  return units_impl::quantity_t::from_value(a);				\
     }

#else

/*!
 * Defines two user-defined literals operators to use with units.
 *
 * \param suffix The suffix you want to write in the code
 * \param qty_type Type of variable (must be a boost::unit::quantity)xk

 */
#define MAKE_LITERAL_OPERATOR(suffix, qty_type)				\
     inline units_impl::qty_type operator"" suffix(const long double a)	\
     {									\
	  return units_impl::qty_type::from_value(a);			\
     }									\
     inline units_impl::qty_type operator"" suffix(const unsigned long long a) \
     {									\
	  return units_impl::qty_type::from_value(a);			\
     }

/*!
 * Defines two user-defined literals operators to use with Boost units.
 *
 * \param suffix The suffix you want to write in the code
 * \param type The type of the corresponding unit.
 * \param value Boost Units value of the \c type unit
 */
#define MAKE_LITERAL_OPERATOR_BOOST(suffix, type, value)		\
     inline boost::units::quantity<type>				\
     operator"" suffix(const long double a)				\
     {									\
	  return boost::units::quantity<type> (a * value);		\
     }									\
     inline boost::units::quantity<type>				\
     operator"" suffix(const unsigned long long a)			\
     {									\
	  return boost::units::quantity<type> (a * value);		\
     }

#endif // NO_UNITS

#endif //UNITS_IMPL_HPP_INCLUDED
