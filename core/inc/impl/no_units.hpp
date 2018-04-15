/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef NO_UNITS_HPP_INCLUDED
#define NO_UNITS_HPP_INCLUDED

#include "impl/units_impl.hpp"

#include <ostream>

#include <boost/operators.hpp>

namespace cellsim {
namespace units_impl {
     class quantity_t : 
	boost::arithmetic<quantity_t,
	boost::less_than_comparable<quantity_t,
	boost::equality_comparable<quantity_t,
	boost::unit_steppable<quantity_t,
	boost::arithmetic<quantity_t,double
			  >>>>>
     {
	  typedef double value_type;
     public:
	  constexpr quantity_t() : v_(0.) {}
	  constexpr quantity_t(value_type v) : v_(v) {}


	  value_type value() const {return v_;}

	  static quantity_t from_value(value_type v) {return quantity_t(v);}

	  bool operator==(const quantity_t& a) const
	       { return v_ == a.v_; }

	  quantity_t operator-() const
	  { 
	       quantity_t tmp(*this);
	       tmp.v_ = -tmp.v_;
	       return tmp;
	  }
	  
	  quantity_t& operator++()
	       { ++v_; return *this; }
	  quantity_t& operator--()
	       { --v_; return *this; }

	  bool operator< (const quantity_t& a) const
	       { return v_ < a.v_; }

	  quantity_t& operator+=(const quantity_t& a)
	       { 
		    v_ += a.v_;
		    return *this;
	       }
	  quantity_t& operator-=(const quantity_t& a)
	       {
		    v_ -= a.v_;
		    return *this;
	       }
	  quantity_t& operator*=(const quantity_t& a)
	       {
		    v_ *= a.v_;
		    return *this;
	       }
	  quantity_t& operator/=(const quantity_t& a)
	       {
		    v_ /= a.v_;
		    return *this;
	       }

	  explicit operator value_type() { return v_; }

	  quantity_t& operator+=(double a)
	       { 
		    v_ += a;
		    return *this;
	       }
	  quantity_t& operator-=(double a)
	       {
		    v_ -= a;
		    return *this;
	       }
	  quantity_t& operator*=(double a)
	       {
		    v_ *= a;
		    return *this;
	       }
	  quantity_t& operator/=(double a)
	       {
		    v_ /= a;
		    return *this;
	       }

     private:
	  value_type v_;
     };

     typedef quantity_t concentration_ut;
     typedef quantity_t electric_potential_ut;
     typedef quantity_t potential_SR_ut;
     typedef quantity_t length_ut;
     typedef quantity_t flux_ut;
     typedef quantity_t time_ut;
     typedef quantity_t inverse_time_ut;
     typedef quantity_t diff_coef_ut;
     typedef quantity_t permeablity_ut;

     typedef quantity_t conductance_ut;
     typedef quantity_t gamma_ut;
     typedef quantity_t beta_ut;
} // namespace units_impl

     MAKE_LITERAL_OPERATOR(_s)
     MAKE_LITERAL_OPERATOR(_um)
     MAKE_LITERAL_OPERATOR(_M)
     MAKE_LITERAL_OPERATOR(_uM)
     MAKE_LITERAL_OPERATOR(_mV)

     template <std::size_t power>
     inline units_impl::quantity_t pow(units_impl::quantity_t q)
     {
	  return q * pow<power-1>(q);
     }
     template <>
     inline units_impl::quantity_t pow<0>(units_impl::quantity_t )
     {
	  return 1.;
     }
} // namespace cellsim

     inline std::ostream& operator<<(std::ostream& out, cellsim::units_impl::quantity_t q)
     {
	  out << q.value();
	  return out;
     }

     inline cellsim::units_impl::quantity_t exp(cellsim::units_impl::quantity_t q) {return exp(q.value());}
     

namespace boost {
namespace units {
namespace si {
     struct time;
     constexpr cellsim::units_impl::quantity_t um;
     constexpr cellsim::units_impl::quantity_t uM;
     constexpr cellsim::units_impl::quantity_t mV;
     constexpr cellsim::units_impl::quantity_t second;
     constexpr cellsim::units_impl::quantity_t seconds;
} // namespace si
} // namespace units
} // namespace boost



#endif //NO_UNITS_HPP_INCLUDED
