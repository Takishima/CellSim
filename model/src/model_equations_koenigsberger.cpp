#include "model_equations_koenigsberger.hpp"
#include "cell_koenigsberger.hpp"

namespace model_equations = cellsim::model_koenigsberger::model_equations;
using namespace cellsim::model_koenigsberger::model_functions;

cellsim::flux_ut model_equations::dc_dt(const Cell cell, 
					Constants_f_param_t C,
					cellsim::time_ut /*t*/)
{
     return JIP3(cell.I(), C) 
	  - JSRuptake(cell.c(), C) 
	  + JCICR(cell.c(), cell.s(), C)
	  - Jextrusion(cell.c(), cell.v(), C)
	  + Jleak(cell.s(), C) 
	  - JVOCC(cell.v(), C) 
	  + JNaCa(cell.c(), cell.v(), C);
}

cellsim::flux_ut model_equations::dI_dt(const Cell cell, 
					Constants_f_param_t C,
					cellsim::time_ut /*t*/)
{
     //using constants::JPLCago_bg;

     return JPLCago()
	  - Jdegrad(cell.I(), C);
     // return JPLCago_bg
     // 	  - Jdegrad(cell.I());
}

cellsim::flux_ut model_equations::ds_dt(const Cell cell, 
					Constants_f_param_t C,
					cellsim::time_ut /*t*/)
{
     return JSRuptake(cell.c(), C)
	  - JCICR(cell.c(), cell.s(), C) 
	  - Jleak(cell.s(), C);
}

cellsim::potential_SR_ut model_equations::dv_dt(const Cell cell, 
						Constants_f_param_t C,
						cellsim::time_ut /*t*/)
{
     return C.gamma * (-JNaK(C)
		     - JCl(cell.c(), cell.v(), C)
		     - 2. * JVOCC(cell.v(), C)
		     - JNaCa(cell.c(), cell.v(), C)
		       - JK(cell.v(), cell.w(), C)
		     - Jback(cell.v(), C))
	  + cell.Vcoupling();
	  // + Vcoupling(cell, C);
}

cellsim::inverse_time_ut model_equations::dw_dt(const Cell cell,
						Constants_f_param_t C,
						cellsim::time_ut /*t*/)
{
     using boost::units::si::second;
     return C.lambda * (Kactiv(cell.c(), cell.v(), C) - cell.w());
}
