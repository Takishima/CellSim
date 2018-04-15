#include "model_functions_2009.hpp"

#include <boost/units/pow.hpp>

namespace model_functions = cellsim::model_halidi::model_functions_2009;

using cellsim::flux_ut;
using cellsim::potential_SR_ut;

flux_ut model_functions::JVOCC(electric_potential_ut v, 
			       Constants_f_param_t C) 
{
     return C.GCa * (v - C.vCa1) / (1 + exp(-(v - C.vCa2) / C.RCa));
}

flux_ut model_functions::JNaCa(concentration_ut c, 
			       electric_potential_ut v,
			       Constants_f_param_t C)
{
     return C.GNaCa * c / (c + C.cNaCa) * (v - C.vNaCa);
}

flux_ut model_functions::JSRuptake(concentration_ut c, 
				   Constants_f_param_t C)
{
     return C.B * c * c / (c * c + C.cb * C.cb);
}

flux_ut model_functions::JCICR(concentration_ut c, 
			       concentration_ut s, 
			       Constants_f_param_t C)
{
     using boost::units::pow;
     return C.C * s * s / (C.sc *C.sc + s * s) 
	  * ( pow<4>(c) / (pow<4>(C.cc) + pow<4>(c)) );
}

flux_ut model_functions::Jextrusion(concentration_ut c, 
				    electric_potential_ut v, 
				    Constants_f_param_t C)
{
     return C.D * c * (1 + (v - C.vd) / C.Rd);
}

flux_ut model_functions::Jleak(concentration_ut s, Constants_f_param_t C)
{
     return C.L * s;
}

flux_ut model_functions::JPLCd(concentration_ut c, Constants_f_param_t C)
{
     return C.E * c * c / (C.KCa * C.KCa + c * c);
}

flux_ut model_functions::Jdegrad(concentration_ut I, Constants_f_param_t C)
{
     return C.k * I;
}

flux_ut model_functions::JIP3(concentration_ut I, Constants_f_param_t C)
{
     return C.F * I * I / (C.KI * C.KI + I * I);
}

flux_ut model_functions::JNaK(Constants_f_param_t C)
{
     return C.FNaK;
}

flux_ut model_functions::JCl(electric_potential_ut v,
			     Constants_f_param_t C)
{
     return C.GCl * (v - C.vCl);
}

flux_ut model_functions::JK(electric_potential_ut v, 
			    double w, 
			    Constants_f_param_t C)
{
     return C.GK * w * (v - C.vK);
}

double model_functions::Kactiv(concentration_ut c,
			       electric_potential_ut v, 
			       Constants_f_param_t C)
{
     using boost::units::pow;
     return double(pow<2>(c+C.cw)/(pow<2>(c+C.cw) + C.beta *exp(-(v - C.vCa3) / C.RK)));
}

