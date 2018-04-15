#include "cell.hpp"
#include "output_streams.hpp"
#include "model_functions_2009.hpp"
#include "model_functions_2012.hpp"

#include <ostream>
#include <sstream>

namespace halidi = cellsim::model_halidi;
using boost::units::si::seconds;

halidi::Cell::Cell(size_type n_point, Constants_f_param_t constants)
     : Cell(n_point, 
	    concentration_t(), 
	    concentration_t(), 
	    concentration_t(), 
	    electric_potential_ut(), 
	    0.,
	    constants)
{}

halidi::Cell::Cell(size_type n_point, 
		   concentration_t c_in, 
		   concentration_t I_in,
		   concentration_t s_in,
		   electric_potential_ut v_in,
		   double w_in,
		   Constants_f_param_t constants)
     : n_point_(n_point),
       cell_points_(n_point, c_in, I_in, s_in),
       v_now_(v_in), 
       w_now_(w_in), 
       v_new_(v_in),
       w_new_(w_in),
       prev_(nullptr), next_(nullptr),
       model_constants_(constants)
{
     /*
      * We actually store n_point+2 points inside the cell :
      * - n_point actual points
      * - first point corresponds to uj-1*
      * - last point corresponds to ui+1*
      * (see compute_gap_junctions())
      */
     assert(w_now_ > 0.);

     /*
      * This was used by Salathe in 2009 to get access to the membrane potential
      * of neighbouring cells without storing pointers to them.
      */
     // cell_points_.set_s_bc(concentration_ut::from_value(v_now_.value()), 
     // 			   concentration_ut::from_value(v_now_.value()));
}


// DEPRECATED constructor
halidi::Cell::Cell(size_type n_point, 
		   concentration_t c_in,
		   concentration_t I_in,
		   concentration_t s_in,
		   electric_potential_ut v_in, 
		   double w_in, 
		   permeablity_ut /*Pip_in*/,
		   permeablity_ut /*Pc_in*/,
		   diff_coef_ut /*Dip_in*/,
		   diff_coef_ut /*Dc_in*/,
		   Constants_f_param_t constants)
     : n_point_(n_point),
       cell_points_(n_point, c_in, I_in, s_in),
       v_now_(v_in), 
       w_now_(w_in), 
       v_new_(v_in),
       w_new_(w_in),
       prev_(nullptr), next_(nullptr),
       model_constants_(constants)
{}

// ==============================================================================

void halidi::Cell::set_neighbours(const Cell* prev, const Cell* next)
{
     prev_ = prev;
     next_ = next;
}

/*
time_integrator<V_EQN>::do_step(v, v_, t, dt,
-                                         std::bind(&mod_eq::dv_dt,
-                                                   c_avg(),
-                                                   _1,
-                                                   w,
-                                                   cell_points_.s_bc_front_now() = vm
-                                                   cell_points_.s_bc_back_now() = vp
-                                                   _2),
-                                         activator);


*/

cellsim::potential_SR_ut halidi::Cell::Vcoupling() const
{
     potential_SR_ut res;

     // res += -C.g * (v_now_ - electric_potential_ut::from_value(cell_points_.s_bc_front_now().value()));
     // res += -C.g * (v_now_ - electric_potential_ut::from_value(cell_points_.s_bc_back_now().value()));
     
     if (prev_ != nullptr) {
     	  res += -model_constants_.g * (v_now_ - prev_->v_now_);
     }
     if (next_ != nullptr) {
     	  res += -model_constants_.g * (v_now_ - next_->v_now_);
     }
     
     return res;
}

void halidi::Cell::compute_gap_junctions(length_ut dx)
{
     assert(dx.value() > 0.);
     
     /*
      * To understand what is done here, you should read the report by
      * Nguyen Damien, EPFL, 2012 or have a look at the end of
      * J. Sneyd, B. T. Wetton, A. C. Charles, and M. J. Sanderson,
      *   Intercellular calcium waves mediated by diffusion of inositol 
      *   trisphosphate: a two-dimensional model,
      *   American Journal of Physiology - Cell Physiology 268 (1995)
      */
     /*
      * Note:
      * - first point corresponds to uj-1*
      * - last point corresponds to uj+1*
      */

     /* 
      * Calcium and IP3 diffusion
      * See equations A.IV and A.V in Nguyen's report.
      */
     const auto beta_c(model_constants_.Pc * dx /
		       (model_constants_.Dc + model_constants_.Pc * dx));
     const auto beta_I(model_constants_.Pip * dx / 
		       (model_constants_.Dip + model_constants_.Pip * dx));

     // Calculate uj-1*
     if (prev_ != nullptr) {
	  cell_points_.ca_bc_front_next() = 
	       beta_c * prev_->cell_points_.ca_back_now()
	       + (1 - beta_c) * cell_points_.ca_front_now();

	  cell_points_.ip3_bc_front_next() = 
	       beta_I * prev_->cell_points_.ip3_back_now()
	       + (1 - beta_I) * cell_points_.ip3_front_now();

	  // cell_points_.s_bc_front_next() = concentration_ut::from_value(prev_->v_now_.value());
     }
     else {     
	  cell_points_.ca_bc_front_next() = 
	       beta_c * cell_points_.ca_front_now()
	       + (1 - beta_c) * cell_points_.ca_front_now();
	  
	  cell_points_.ip3_bc_front_next() = 
	       beta_I * cell_points_.ip3_front_now()
	       + (1 - beta_I) * cell_points_.ip3_front_now();

	  // cell_points_.s_bc_front_next() = concentration_ut::from_value(v_now_.value());
     }
     cell_points_.ca_bc_front_now() = 
     	  cell_points_.ca_bc_front_next();
     cell_points_.ip3_bc_front_now() = 
     	  cell_points_.ip3_bc_front_next();
     // cell_points_.s_bc_front_now() = 
     // 	  cell_points_.s_bc_front_next();

     // Calculate ui+1*
     if (next_ != nullptr) {
	  cell_points_.ca_bc_back_next() = 
	       beta_c * next_->cell_points_.ca_front_now()
	       + (1 - beta_c) * cell_points_.ca_back_now();

	  cell_points_.ip3_bc_back_next() =
	       beta_I * next_->cell_points_.ip3_front_now()
	       + (1 - beta_I) * cell_points_.ip3_back_now();

	  // cell_points_.s_bc_back_next() = concentration_ut::from_value(next_->v_now_.value());
     }
     else {
	  cell_points_.ca_bc_back_next() = 
	       beta_c * cell_points_.ca_back_now()
	       + (1 - beta_c) * cell_points_.ca_back_now();

	  cell_points_.ip3_bc_back_next() =
	       beta_I * cell_points_.ip3_back_now()
	       + (1 - beta_I) * cell_points_.ip3_back_now();

	  // cell_points_.s_bc_back_next() = concentration_ut::from_value(v_now_.value());
     }
     cell_points_.ca_bc_back_now() = 
     	  cell_points_.ca_bc_back_next();
     cell_points_.ip3_bc_back_now() = 
     	  cell_points_.ip3_bc_back_next();
     // cell_points_.s_bc_back_now() = 
     // 	  cell_points_.s_bc_back_next();
}

void halidi::Cell::swap_times()
{
     cell_points_.swap_times();
     std::swap(v_now_, v_new_);
     std::swap(w_now_, w_new_);
}

void halidi::Cell::print_cell(output::Output_Streams& out,
			      time_ut t,
			      size_type& point_offset) const
{
     const auto END(cend<CA, NOW>(cell_points_));
     for(auto c_it(cbegin<CA, NOW>(cell_points_)),
	      I_it(cbegin<IP3, NOW>(cell_points_)),
	      s_it(cbegin<S, NOW>(cell_points_)) 
	      ; c_it != END
	      ; ++c_it, ++I_it, ++s_it) {

	  out.ca_out << t.value() << " "
		     << point_offset << " " 
		     << c_it->value() << std::endl;

	  out.ip3_out << t.value() << " "
		      << point_offset << " " 
		      << I_it->value() << std::endl;

	  out.s_out << t.value() << " "
		    << point_offset << " " 
		    << s_it->value() << std::endl;

	  ++point_offset;
     }
     out.v_out << v_now_.value() << std::endl;
     out.w_out << w_now_ << std::endl;
}

void halidi::Cell::print_point(output::Output_Streams& out,
			       size_type point_idx) const
{
     out.ca_out << cell_points_.get<CA>(point_idx).value() << " ";
     out.ip3_out << cell_points_.get<IP3>(point_idx).value() << " ";
     out.s_out << cell_points_.get<S>(point_idx).value() << " ";
     out.v_out << v_now_.value() << " ";
     out.w_out << w_now_ << " ";
}

std::string halidi::Cell::get_fluxes_header()
{
     std::string one("# t JVOCC JNaCa JNaK JSRuptake JCICR Jextrusion Jleak Jdegrad ");
     std::string two("# 1   2     3     4     5        6       7        8      9    ");

#ifdef HALIDI_2009
     one += "JIP3 JCl ";
     two += " 10   11 ";
#else
     one += "JIP3 JCl Jback ";
     two += " 10   11   12  ";
#endif

     one += "JK Kactiv Vcoupling\n";
#ifdef HALIDI_2009
     two += "12   13       14   \n";
#else
     two += "13   14       15   \n";
#endif
     return (one + two);
}

void halidi::Cell::print_fluxes(std::ostream& out,
				time_ut t,
				size_type point_idx) const
{
#ifdef HALIDI_2009
     using namespace model_functions_2009;
#else
     using namespace model_functions_2012;
#endif // HALIDI_2009

     typedef electric_potential_ut e_p_ut;

     const auto& c = cell_points_.get<CA>(point_idx);
     const auto& I = cell_points_.get<IP3>(point_idx);
     const auto& s = cell_points_.get<S>(point_idx);
     
     out << t.value() << " "
	 << JVOCC(v_now_, model_constants_).value() << " "
	 << JNaCa(c, v_now_, model_constants_).value() << " "
	 << JNaK(model_constants_).value() << " "
	 << JSRuptake(c, model_constants_).value() << " "
	 << JCICR(c, s, model_constants_).value() << " "
	 << Jextrusion(c, v_now_, model_constants_).value() << " "
	 << Jleak(s, model_constants_).value() << " "
	 << Jdegrad(I, model_constants_).value() << " "
	 << JIP3(I, model_constants_).value() << " "
#ifdef HALIDI_2009
	 << JCl(v_now_, model_constants_).value() << " "
#else
	 << JCl(c, v_now_, model_constants_).value() << " "
	 << Jback(v_now_, model_constants_).value() << " "
#endif
	 << JK(v_now_, w_now_, model_constants_).value() << " "
	 << Kactiv(c, v_now_, model_constants_) << " "
	 << Vcoupling().value() << " ";
}

#ifdef TEST
namespace ba = boost::accumulators;
typedef ba::accumulator_set<double, ba::stats<ba::tag::mean, 
					      ba::tag::min,
					      ba::tag::max,
					      ba::tag::variance> > acc_type;

std::tuple<double, double, double, double> halidi::Cell::getca() const 
{
     typedef container_type::value_type value_t;
     using namespace boost::accumulators;

     acc_type acc;
     std::for_each(cbegin<CA, NOW>(cell_points_),
		   cend<CA, NOW>(cell_points_),
		   [&] (concentration_ut c) {
			acc(c.value());
		   });

     return std::make_tuple(mean(acc), variance(acc), ba::min(acc), ba::max(acc));
}
	  
std::tuple<double, double, double, double> halidi::Cell::getI() const 
{
     typedef container_type::value_type value_t;
     using namespace boost::accumulators;

     acc_type acc;
     std::for_each(cbegin<IP3, NOW>(cell_points_),
		   cend<IP3, NOW>(cell_points_),
		   [&] (concentration_ut c) {
			acc(c.value());
		   });

     return std::make_tuple(mean(acc), variance(acc), ba::min(acc), ba::max(acc));
}

std::tuple<double, double, double, double> halidi::Cell::gets() const
{
     typedef container_type::value_type value_t;
     using namespace boost::accumulators;

     acc_type acc;
     std::for_each(cbegin<S, NOW>(cell_points_),
		   cend<S, NOW>(cell_points_),
		   [&] (concentration_ut c) {
			acc(c.value());
		   });

     return std::make_tuple(mean(acc), variance(acc), ba::min(acc), ba::max(acc));
}

std::tuple<double, double, double, double, double, double> halidi::Cell::test_compute_gap_junctions(length_ut dx)
{
     compute_gap_junctions(dx);

     auto ca_back = cell_points_.ca_bc_back_next();
     auto ca_front = cell_points_.ca_bc_front_next();

     auto I_front = cell_points_.ip3_bc_front_next();
     auto I_back = cell_points_.ip3_bc_back_next();

     electric_potential_ut s_front;
     if (prev_ == nullptr) {
	  s_front = v_now_;
     }
     else {
	  s_front = prev_->v_now_;
     }
     electric_potential_ut s_back;
     if (next_ == nullptr) {
	  s_back = v_now_;
     }
     else {
	  s_back = next_->v_now_;
     }

     return std::make_tuple(ca_front.value(), I_front.value(), s_front.value(),
			    ca_back.value(), I_back.value(), s_back.value());
}


cellsim::potential_SR_ut halidi::Cell::Vcoupling(Constants_f_param_t C,
						 electric_potential_ut vm,
						 electric_potential_ut vp) const
{
     potential_SR_ut res;

     res += -C.g * (v_now_ - vm);
     res += -C.g * (v_now_ - vp);
     
     // if (prev_ != nullptr) {
     // 	  res += -C.g * (v_now_ - prev_->v_now_);
     // }
     // if (next_ != nullptr) {
     // 	  res += -C.g * (v_now_ - next_->v_now_);
     // }
     
     return res;
}

#endif // TEST
