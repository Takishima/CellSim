#include "model_functions_2012.hpp"

namespace functions = cellsim::model_halidi::model_functions_2012;

using cellsim::flux_ut;
using cellsim::potential_SR_ut;


flux_ut functions::JVOCC(electric_potential_ut v, Constants_f_param_t C)
{
     return C.GCa * (v - C.vCa1) / (1 + exp(-(v - C.vCa2) / C.RCa));
}

flux_ut functions::JNaCa(concentration_ut c, 
			 electric_potential_ut v, 
			 Constants_f_param_t C) 
{
     return C.GNaCa * (c / (c + C.cNaCa)) * (v - C.vNaCa);
}

flux_ut functions::JSRuptake(concentration_ut c, Constants_f_param_t C) 
{
     return C.B * c * c / (c * c + C.cb * C.cb);
}

flux_ut functions::JCICR(concentration_ut c,
			 concentration_ut s,
			 Constants_f_param_t C) 
{
     const auto c4(c * c * c * c);
     return C.C * (s * s / (C.sc * C.sc + s * s))
	  * (c4 / (C.cc * C.cc * C.cc * C.cc + c4));
}

flux_ut functions::Jextrusion(concentration_ut c, 
			      electric_potential_ut v, 
			      Constants_f_param_t C)
{
     return C.D * c * (1 + (v - C.vd) / C.Rd);
}

flux_ut functions::Jleak(concentration_ut s, Constants_f_param_t C) 
{
     return C.L * s;
}

flux_ut functions::JPLCd(concentration_ut c, Constants_f_param_t C) 
{
     return C.E * c * c / (C.KCa * C.KCa + c * c);
}

flux_ut functions::Jdegrad(concentration_ut I, Constants_f_param_t C) 
{
     return C.k * I;
}

flux_ut functions::JIP3(concentration_ut I, Constants_f_param_t C) 
{
     return C.F * I * I / (C.KI * C.KI + I * I);
}

flux_ut functions::JNaK(Constants_f_param_t C) 
{
     return C.FNaK;
}

flux_ut functions::JCl(concentration_ut c, 
		       electric_potential_ut v, 
		       Constants_f_param_t C)
{
     return C.GCl * (c / (c + C.cCl)) * (v - C.vCl);
}

flux_ut functions::Jback(electric_potential_ut v, Constants_f_param_t C)
{
     return C.Gback * (v - C.vrest);
}

flux_ut functions::JK(electric_potential_ut v, double w, Constants_f_param_t C) 
{
     return C.GK * w * (v - C.vK);
}

double functions::Kactiv(concentration_ut c,
			 electric_potential_ut v, 
			 Constants_f_param_t C)
{
     const auto sum(c + C.cw);
     return double(sum * sum / (sum * sum + C.beta * exp(-(v - C.vCa3) / C.RK)));
}

