#include <algorithm>
#include <cassert>
#include <cmath>

#include "model_functions_koenigsberger.hpp"
#include "cell_koenigsberger.hpp"
#include "constants_koenigsberger.hpp"

namespace functions = cellsim::model_koenigsberger::model_functions;

cellsim::flux_ut functions::JVOCC(electric_potential_ut v, 
				  Constants_f_param_t C)
{
     return C.GCa * (v - C.vCa1) / (1 + exp(-(v - C.vCa2) / C.RCa));
}

cellsim::flux_ut functions::JNaCa(concentration_ut c, 
				  electric_potential_ut v, 
				  Constants_f_param_t C) 
{
     return C.GNaCa * (c / (c + C.cNaCa)) * (v - C.vNaCa);
}

cellsim::flux_ut functions::JSRuptake(concentration_ut c, 
				      Constants_f_param_t C) 
{
     return C.B * c * c / (c * c + C.cb * C.cb);
}

cellsim::flux_ut functions::JCICR(concentration_ut c, 
				  concentration_ut s, 
				  Constants_f_param_t C) 
{
     const auto c4(c * c * c * c);

     return C.C * (s * s / (C.sc * C.sc + s * s))
	  * (c4 / (C.cc * C.cc * C.cc * C.cc + c4));
}

cellsim::flux_ut functions::Jextrusion(concentration_ut c,
				       electric_potential_ut v, 
				       Constants_f_param_t C)
{
     return C.D * c * (1 + (v - C.vd) / C.Rd);
}

cellsim::flux_ut functions::Jleak(concentration_ut s, 
				  Constants_f_param_t C) 
{
     return C.L * s;
}

cellsim::flux_ut functions::Jdegrad(concentration_ut I, 
				  Constants_f_param_t C) 
{
     return C.k * I;
}

cellsim::flux_ut functions::JIP3(concentration_ut I, 
				 Constants_f_param_t C) 
{
     return C.F * I * I / (C.KI * C.KI + I * I);
}

cellsim::flux_ut functions::JNaK(Constants_f_param_t C) 
{
     return C.FNaK;
}

cellsim::flux_ut functions::JCl(concentration_ut c, 
				electric_potential_ut v, 
				Constants_f_param_t C)
{
     return C.GCl * (c / (c + C.cCl)) * (v - C.vCl);
}

cellsim::flux_ut functions::Jback(electric_potential_ut v, 
				  Constants_f_param_t C)
{
     return C.Gback * (v - C.vrest);
}

cellsim::flux_ut functions::JK(electric_potential_ut v, 
			       double w, 
			       Constants_f_param_t C) 
{
     return C.GK * w * (v - C.vK);
}

double functions::Kactiv(concentration_ut c, 
			 electric_potential_ut v, 
			 Constants_f_param_t C)
{
     concentration_ut sum(c + C.cw);
     return double(sum * sum / (sum * sum + C.beta * exp(-(v - C.vCa3) / C.RK)));
}

cellsim::potential_SR_ut functions::Vcoupling(const Cell c, 
					      Constants_f_param_t C)
{
     electric_potential_ut result(0._mV);
     for(const auto other : c.neighbours()) {
	  result += (c.v() - other->v());
     }
     return -C.g * result;
}

