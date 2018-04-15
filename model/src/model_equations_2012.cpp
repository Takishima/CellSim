#include "model_equations_2012.hpp"
#include "cell.hpp"

namespace model_equations = cellsim::model_halidi::model_equations_2012;
using namespace cellsim::model_halidi::model_functions_2012;
using cellsim::flux_ut;
using cellsim::potential_SR_ut;
using cellsim::inverse_time_ut;

flux_ut model_equations::dc_dt(const Point_Index p, 
			       electric_potential_ut v,
			       Constants_f_param_t C,
			       time_ut /*t*/)
{
     return JIP3(p.get<IP3>(), C) 
	  - JSRuptake(p.get<CA>(), C) 
	  + JCICR(p.get<CA>(), p.get<S>(), C)
	  - Jextrusion(p.get<CA>(), v, C)
	  + Jleak(p.get<S>(), C)
	  - JVOCC(v, C) 
	  + JNaCa(p.get<CA>(), v, C);
}

flux_ut model_equations::dI_dt(const Point_Index p, 
			       flux_ut JPLCago,
			       Constants_f_param_t C,
			       time_ut /*t*/)
{
     return JPLCago
	  + JPLCd(p.get<CA>(), C) 
	  - Jdegrad(p.get<IP3>(), C);
}

flux_ut model_equations::ds_dt(const Point_Index p,
			       Constants_f_param_t C,
			       time_ut /*t*/)
{
     return JSRuptake(p.get<CA>(), C)
	  - JCICR(p.get<CA>(), p.get<S>(), C) 
	  - Jleak(p.get<S>(), C);
}

potential_SR_ut model_equations::dv_dt(const Cell c,
				       Constants_f_param_t C,
				       time_ut /*t*/)
{
     return C.gamma * (-JNaK(C)
		       - JCl(c.c_avg(), c.v(), C)
		       - 2. * JVOCC(c.v(), C)
		       -JNaCa(c.c_avg(), c.v(), C)
		       - JK(c.v(), c.w(), C)
		       - Jback(c.v(), C))
	  + c.Vcoupling();
}

inverse_time_ut model_equations::dw_dt(const Cell c,
				       Constants_f_param_t C,
				       time_ut /*t*/)
{
     return C.lambda * (Kactiv(c.c_avg(), c.v(), C) - c.w());
}
