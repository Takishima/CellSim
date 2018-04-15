#include "constants_koenigsberger.hpp"

namespace model = cellsim::model_koenigsberger;

cellsim::flux_ut& model::JPLCago()
{
     static flux_ut jplcago_(model::constants::JPLCago_bg);
     return jplcago_;
}
