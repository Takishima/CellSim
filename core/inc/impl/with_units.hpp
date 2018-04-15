/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef WITH_UNITS_HPP_INCLUDED
#define WITH_UNITS_HPP_INCLUDED

#include "impl/units_impl.hpp"

#include <boost/units/io.hpp>
#include <boost/units/dimension.hpp>
#include <boost/units/make_scaled_unit.hpp>
#include <boost/units/reduce_unit.hpp>
  
#include <boost/units/systems/si/codata_constants.hpp>
#include <boost/units/systems/si/amount.hpp>
#include <boost/units/systems/si/electric_potential.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/prefixes.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/volume.hpp>

#pragma GCC system_header

namespace cellsim {
//! Units implementation namespace
namespace units_impl {
     namespace mpl = boost::mpl;
     namespace bu = boost::units;
     namespace si = bu::si;

     typedef mpl::divides<bu::amount_dimension, 
			  bu::volume_dimension>::type concentration_dimension;
     typedef mpl::divides<concentration_dimension,
			  bu::time_dimension>::type flux_dimension;

     // Define 1 molar = 1 mol / L = 1kmol / m3
     typedef bu::make_scaled_unit<
	  bu::unit<concentration_dimension, si::system>, 
	  bu::scale<10, bu::static_rational<3> > >::type molar_unit;

     // Define 1 micro molar
     typedef bu::make_scaled_unit<
	  molar_unit, 
	  bu::scale<10, bu::static_rational<-6> > >::type micro_molar_unit;

     // Define 1 milli volt
     typedef bu::make_scaled_unit<
	  si::electric_potential, 
	  bu::scale<10, bu::static_rational<-3> > >::type milli_volt_unit;

     // Define 1 micro meter
     typedef bu::make_scaled_unit<
	  si::length, 
	  bu::scale<10, bu::static_rational<-6> > >::type micro_metre_unit;

     // Define flux units
     typedef bu::unit<flux_dimension, si::system> flux_unit;
     typedef bu::make_scaled_unit<
	  flux_unit,
	  bu::scale<10, bu::static_rational<3> > >::type molar_per_second_unit;
     typedef bu::make_scaled_unit<
	  molar_per_second_unit, 
	  bu::scale<10, bu::static_rational<-6> > >::type micro_molar_per_second_unit;

     // Define potential slew rate unit (similar to V/s)
     typedef mpl::divides<bu::electric_potential_dimension, 
			  bu::time_dimension>::type potential_slew_rate;
     typedef bu::unit<potential_slew_rate, si::system> potential_SR_unit_base;
     typedef bu::make_scaled_unit<
	  potential_SR_unit_base, 
	  bu::scale<10, bu::static_rational<-3> > >::type potential_SR_unit;


     // Define diffusion coefficient dimension in [um^2 / s]
     typedef mpl::divides<
	  mpl::multiplies<bu::length_dimension, 
			  bu::length_dimension>::type,
	  bu::time_dimension>::type diffusion_dimension;
     typedef bu::unit<diffusion_dimension, si::system> diffusion_unit_base;
     typedef bu::make_scaled_unit<
	  diffusion_unit_base, 
	  bu::scale<10, bu::static_rational<-12> > >::type diffusion_unit;

     // Define permeability dimension in [um / s]
     typedef mpl::divides<bu::length_dimension, 
			  bu::time_dimension>::type permeablity_dimension;
     typedef bu::unit<permeablity_dimension, si::system> permeability_unit_base;
     typedef bu::make_scaled_unit<
	  permeability_unit_base, 
	  bu::scale<10, bu::static_rational<-6> > >::type permeablity_unit;

     // Define inverse time dimension
     typedef bu::derived_dimension<bu::time_base_dimension,-1>::type inverse_time_dimension;

     // Quantity typedefs
     /*!
      * \brief Concentration quantities
      * Units are micro molars.
      */
     typedef bu::quantity<micro_molar_unit> concentration_ut;
     /*!
      * \brief Electric potential quantities
      * Units are milli volts.
      */
     typedef bu::quantity<milli_volt_unit> electric_potential_ut;
     /*!
      * \brief Distance/Length quantities
      * Units are micro metres.
      */
     typedef bu::quantity<micro_metre_unit> length_ut;
     /*!
      * \brief Electric potential Slew Rate quantities
      * Units are milli volts per seconds
      */
     typedef bu::quantity<potential_SR_unit> potential_SR_ut;
     /*!
      * \brief Flux quantities
      * Units are micro molar per seconds.
      */
     typedef bu::quantity<micro_molar_per_second_unit> flux_ut;
     /*!
      * \brief Time quantities
      * Units are seconds
      */
     typedef bu::quantity<boost::units::si::time> time_ut;
     /*!
      * \brief Inverse time quantities
      * Units are inverse seconds.
      */
     typedef bu::quantity<bu::unit<inverse_time_dimension, 
				   si::system>> inverse_time_ut;
     /*!
      * \brief Diffusion coefficient quantities
      * Units are (micro metre)^2 / s
      */
     typedef bu::quantity<diffusion_unit> diff_coef_ut;
     /*!
      * \brief Permeability coefficient quantities
      * Units are micro metre / s
      */
     typedef bu::quantity<permeablity_unit> permeablity_ut;


} // namespace units_impl

     BOOST_UNITS_STATIC_CONSTANT(s_, boost::units::si::time);

     BOOST_UNITS_STATIC_CONSTANT(molar, units_impl::molar_unit);

     BOOST_UNITS_STATIC_CONSTANT(micro_molar, units_impl::micro_molar_unit);
     BOOST_UNITS_STATIC_CONSTANT(uM, units_impl::micro_molar_unit);

     BOOST_UNITS_STATIC_CONSTANT(micro_metre, units_impl::micro_metre_unit);
     BOOST_UNITS_STATIC_CONSTANT(um, units_impl::micro_metre_unit);

     BOOST_UNITS_STATIC_CONSTANT(milli_volt, units_impl::milli_volt_unit);
     BOOST_UNITS_STATIC_CONSTANT(mV, units_impl::milli_volt_unit);

     BOOST_UNITS_STATIC_CONSTANT(micro_molar_per_second, 
				 units_impl::micro_molar_per_second_unit);

     // ====================

     MAKE_LITERAL_OPERATOR_BOOST(_s, 
     				 boost::units::si::time,
     				 boost::units::si::second)
     MAKE_LITERAL_OPERATOR(_um, length_ut)
     // MAKE_LITERAL_OPERATOR(_M, ###)
     MAKE_LITERAL_OPERATOR(_uM, concentration_ut)
     MAKE_LITERAL_OPERATOR(_mV, electric_potential_ut)

     // ====================

namespace units_impl {
     using boost::units::si::seconds;

     // Some remaining definitions (mainly for constants)
     typedef decltype(0.00129 * uM / mV / seconds) conductance_ut;
     typedef decltype(1970 * mV / uM) gamma_ut;
     typedef decltype(0.13 * uM * uM) beta_ut;
} // namespace units_impl
} // namespace cellsim

#endif //WITH_UNITS_HPP_INCLUDED
