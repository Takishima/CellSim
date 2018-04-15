#include "model_equations_2009.hpp"
#include "cell.hpp"

namespace model_equations = cellsim::model_halidi::model_equations_2009;
using namespace cellsim::model_halidi::model_functions_2009;
using cellsim::flux_ut;
using cellsim::potential_SR_ut;
using cellsim::inverse_time_ut;

flux_ut model_equations::dc_dt(const Point_Index p, 
			       electric_potential_ut v, 
			       Constants_f_param_t C,
			       time_ut /*t*/) 
{
     return JCICR(p.get<CA>(), p.get<S>(), C) 
	  - JSRuptake(p.get<CA>(), C)
	  + Jleak(p.get<S>(), C)
	  - Jextrusion(p.get<CA>(),v, C)
	  + JIP3(p.get<IP3>(), C)
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
		       - JCl(c.v(), C) 
		       - 2. * JVOCC(c.v(), C)
		       - JNaCa(c.c_avg(), c.v(), C)
		       - JK(c.v(), c.w(), C))
	  + c.Vcoupling();
	  // + Vcoupling(c.v, vm, vp, C);
}

inverse_time_ut model_equations::dw_dt(const Cell c,
				       Constants_f_param_t C,
				       time_ut /*t*/)
{
     return C.lambda * (Kactiv(c.c_avg(), c.v(), C) - c.w());
}



#ifdef TEST
potential_SR_ut model_equations::dv_dt_test(const Cell c,
					    electric_potential_ut vm,
					    electric_potential_ut vp,
					    Constants_f_param_t C,
					    time_ut /*t*/)
{
     return C.gamma * (-JNaK(C)
		       - JCl(c.v(), C) 
		       - 2. * JVOCC(c.v(), C)
		       - JNaCa(c.c_avg(), c.v(), C)
		       - JK(c.v(), c.w(), C))
	  + c.Vcoupling(C, vm, vp);
	  // + Vcoupling(c.v, vm, vp, C);
}
#endif //TEST


// potential_SR_ut model_equations::dv_dt(concentration_ut c, 
// 				       electric_potential_ut v,
// 				       double w, 
// 				       electric_potential_ut vm,
// 				       electric_potential_ut vp, 
// 				       Constants_f_param_t C,
// 				       time_ut /*t*/)
// {
//      return C.gamma * (-JNaK(C) - JCl(v, C) - 2. * JVOCC(v, C)
// 		     - JNaCa(c, v, C) - JK(v, w, C))
// 	  + Vcoupling(v, vm, vp, C);
// }

// inverse_time_ut model_equations::dw_dt(concentration_ut c, 
// 				       electric_potential_ut v, 
// 				       double w, 
// 				       Constants_f_param_t C,
// 				       time_ut /*t*/)
// {
//      return C.lambda * (Kactiv(c, v, C) - w);
// }

