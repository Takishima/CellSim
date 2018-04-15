#include "cell_koenigsberger.hpp"
#include "output_streams.hpp"

namespace model = cellsim::model_koenigsberger;

model::Cell::Cell(concentration_ut c_in, 
		  concentration_ut I_in,
		  concentration_ut s_in,
		  electric_potential_ut v_in,
		  double w_in,
		  Constants_f_param_t constants)
     : c_now_(c_in), 
       I_now_(I_in), 
       s_now_(s_in), 
       v_now_(v_in), 
       w_now_(w_in), 
       c_new_(c_in),
       I_new_(I_in),
       s_new_(s_in),
       v_new_(v_in),
       w_new_(w_in),
       model_constants_(constants),
       neighbours_()
{}

void model::Cell::set_neighbours(const neighbours_array array)
{
     neighbours_.clear();
     for(auto ptr : array) {
	  if (ptr != nullptr) {
	       neighbours_.push_back(ptr);
	  }
     }
}

cellsim::potential_SR_ut model::Cell::Vcoupling() const
{
     potential_SR_ut res;

     for (auto other : neighbours_) {
	  res += -model_constants_.g * (v_now_ - other->v_now_);
     }

     return res;
}

void model::Cell::print(std::ostream& out, TYPE value) const
{
     if (value == cellsim::CA) {
	  out << c_now_.value() << " ";
     }
     else if (value == cellsim::IP3) {
	  out << I_now_.value() << " ";
     }
     else if (value == cellsim::S) {
	  out << s_now_.value() << " ";
     }
     else if (value == cellsim::V) {
	  out << v_now_.value() << " ";
     }
     else if (value == cellsim::W) {
	  out << w_now_ << " ";
     }
     else {
	  out << "INVALID ";
     }
}

void model::Cell::print(output::Output_Streams& out) const
{
     out.ca_out << c_now_.value() << " ";
     out.ip3_out << I_now_.value() << " ";
     out.s_out << s_now_.value() << " ";
     out.v_out << v_now_.value() << " ";
     out.w_out << w_now_ << " ";
}

std::string model::Cell::get_fluxes_header()
{
     std::string one("# t JVOCC JNaCa JNaK JSRuptake JCICR Jextrusion Jleak Jdegrad ");
     std::string two("# 1   2     3     4     5        6       7        8      9    ");
     one += "JIP3 JCl Jback JK Kactiv Vcoupling\n";
     two += " 10   11   12  13   14       15   \n";

     return (one + two);
}

void model::Cell::print_fluxes(std::ostream& out,
			       time_ut t) const
{
     using namespace model_koenigsberger::model_functions;

     out << t.value() << " "
	 << JVOCC(v_now_, model_constants_).value() << " "
	 << JNaCa(c_now_, v_now_, model_constants_).value() << " "
	 << JNaK(model_constants_).value() << " "
	 << JSRuptake(c_now_, model_constants_).value() << " "
	 << JCICR(c_now_, s_now_, model_constants_).value() << " "
	 << Jextrusion(c_now_, v_now_, model_constants_).value() << " "
	 << Jleak(s_now_, model_constants_).value() << " "
	 << Jdegrad(I_now_, model_constants_).value() << " "
	 << JIP3(I_now_, model_constants_).value() << " "
	 << JCl(c_now_, v_now_, model_constants_).value() << " "
	 << Jback(v_now_, model_constants_).value() << " "
	 << JK(v_now_, w_now_, model_constants_).value() << " "
	 << Kactiv(c_now_, v_now_, model_constants_) << " "
	 << Vcoupling().value() << " ";
}

void model::Cell::swap_times()
{
     std::swap(c_now_, c_new_);
     std::swap(I_now_, I_new_);
     std::swap(s_now_, s_new_);
     std::swap(v_now_, v_new_);
     std::swap(w_now_, w_new_);
}
