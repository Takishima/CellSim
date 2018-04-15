/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef CONSTANTS_KEONIGSBERGER_HPP_INCLUDED
#define CONSTANTS_KEONIGSBERGER_HPP_INCLUDED

#include "config.hpp"
#include "units.hpp"
#include "constants_impl.hpp"

CLANG_DIAG_OFF(global-constructors)
CLANG_DIAG_OFF(unused-variable)

namespace cellsim {
     namespace model_koenigsberger {
	  GCC_DIAG_OFF(padded)
	  GCC_DIAG_OFF(effc++)

	  struct Constants_Values
	  {
	       Constants_Values() = default;
	       Constants_Values(Constants_Values&&) = default;
	       Constants_Values& operator=(Constants_Values&&) = default;
	       Constants_Values(const Constants_Values&) = default;
	       Constants_Values& operator=(const Constants_Values&) = default;
		    
	       flux_ut F;
	       concentration_ut KI;
	       conductance_ut GCa;
	       electric_potential_ut vCa1;
	       electric_potential_ut vCa2;
	       electric_potential_ut RCa;
	       conductance_ut GNaCa;
	       concentration_ut cNaCa;
	       electric_potential_ut vNaCa;
	       flux_ut B;
	       concentration_ut cb;
	       flux_ut C;
	       concentration_ut sc;
	       concentration_ut cc;
	       inverse_time_ut D;
	       electric_potential_ut vd;
	       electric_potential_ut Rd;
	       inverse_time_ut L;
	       gamma_ut gamma;
	       flux_ut FNaK;
	       conductance_ut GCl;
	       concentration_ut cCl;
	       electric_potential_ut vCl;
	       conductance_ut GK;
	       electric_potential_ut vK;
	       inverse_time_ut lambda;
	       concentration_ut cw;
	       beta_ut beta;
	       electric_potential_ut vCa3;
	       electric_potential_ut RK;
	       conductance_ut Gback;
	       electric_potential_ut vrest;
	       inverse_time_ut k;
	       flux_ut E;
	       concentration_ut KCa;
	       inverse_time_ut g;
	  };
	  
	  typedef const Constants_Values& Constants_f_param_t;

	  namespace constants {
	       // template <typename T>
	       // using q = boost::units::quantity<T>;
	       using namespace boost::units::si;

	       STATIC_CONSTANT(F, 1_uM / seconds);
	       // STATIC_CONSTANT(F, 3.45_uM / seconds);
	       STATIC_CONSTANT(Kr, 1._uM);
	       STATIC_CONSTANT(GCa, 0.036195_uM / mV / seconds);
	       STATIC_CONSTANT(vCa1, 100.0_mV);
	       STATIC_CONSTANT(vCa2, -24._mV);
	       STATIC_CONSTANT(RCa, 8.5_mV);
	       STATIC_CONSTANT(GNaCa, 0.006_uM / mV / seconds);
	       STATIC_CONSTANT(cNaCa, 0.5_uM);
	       STATIC_CONSTANT(vNaCa, -30._mV);
	       STATIC_CONSTANT(B, 49.5_uM / seconds);
	       STATIC_CONSTANT(cb, 1._uM);
	       STATIC_CONSTANT(C, 50.0_uM / seconds);
	       // STATIC_CONSTANT(C, 1545.0_uM / seconds);
	       // STATIC_CONSTANT(C, 900.0_uM / seconds);
	       STATIC_CONSTANT(sc, 2._uM);
	       STATIC_CONSTANT(cc, 0.9_uM);
	       STATIC_CONSTANT(D, 3.6 / seconds);
	       STATIC_CONSTANT(vd, -100._mV);
	       STATIC_CONSTANT(Rd, 250._mV);
	       STATIC_CONSTANT(L, 0.375 / seconds);
	       STATIC_CONSTANT(gamma, 492.5_mV / uM);
	       STATIC_CONSTANT(FNaK, 0.03_uM / seconds);
	       STATIC_CONSTANT(GCl, 0.6_uM / mV / seconds);
	       STATIC_CONSTANT(cCl, 0.7_uM);
	       STATIC_CONSTANT(vCl, -25._mV);
	       STATIC_CONSTANT(GK, 0.0045_uM / mV / seconds);
	       STATIC_CONSTANT(vK, -94._mV);
	       STATIC_CONSTANT(lambda, 675.);
	       STATIC_CONSTANT(cw, 0._uM);
	       STATIC_CONSTANT(beta, 0.001_uM * uM);
	       STATIC_CONSTANT(vCa3, -27._mV);
	       STATIC_CONSTANT(RK, 12._mV);
	       STATIC_CONSTANT(Gback, 0.06_uM / mV / seconds);
	       STATIC_CONSTANT(vrest, -55._mV);
	       STATIC_CONSTANT(JPLCago_bg, 0.05_uM / seconds);
	       STATIC_CONSTANT(JPLCago_act, 0.4_uM / seconds);
	       STATIC_CONSTANT(k, 0.1 / seconds);
	       STATIC_CONSTANT(g, 1000 / seconds);
	  } // namespace constants
	  //! Getter/setter on the current value of background  JPLCagonist flux.
	  flux_ut& JPLCago();

	  GCC_DIAG_ON(effc++)
	  GCC_DIAG_ON(padded)
     } // namespace model_koenigsberger
} // namespace cellsim

CLANG_DIAG_ON(unused-variable)
CLANG_DIAG_ON(global-constructors)

#endif //CONSTANTS_KEONIGSBERGER_HPP_INCLUDED
