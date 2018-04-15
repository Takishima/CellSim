#include "system.hpp"
#include "output_streams.hpp"
#include "flux_streams.hpp"
#include "utility.hpp"

namespace halidi = cellsim::model_halidi;

halidi::System::System(size_type Ncells,
		       const Cell c, 
		       length_ut dx_in,
		       bool periodic_bc)
     : t_(0._s), dx_(dx_in), cells_(Ncells, c)
{
     set_neighbours_(periodic_bc);
}

std::string halidi::System::get_system_header()
{
     return std::string("# t idx value\n");
}

void halidi::System::do_step(time_ut dt, flux_ut JPLCago)
{
     using namespace cellsim::equation_type;
     
     integrators::FTCS<CA_EQN> ca_integrator(dt, dx_);
     integrators::FTCS<IP3_EQN> ip3_integrator(dt, dx_);

#ifdef MULTI_THREAD
     size_type start{0}, step{0}, change{0};
     std::tie(step, change) = utility::threads_step(cells_.size(),
						    n_threads());
     std::vector<std::thread> threads;
     const auto END(end(cells_));
     for (auto it(begin(cells_)) ; it != END ; it += step, start += step) {
     	  if (start == change) {
     	       --step;
	  }

#ifndef NDEBUG
	  std::cout << "Launching thread for cells from "
		    << (it - begin(cells_)) << " to "
		    << (it - begin(cells_))+step
		    << std::endl;
#endif // NDEBUG
	  threads.push_back(
	       std::thread([=] ()
			   {
				const auto last(it+step);
				for(auto cell(it); cell != last ; ++cell) {
				     cell->do_step<integrator_t>(t_, dt,
								 JPLCago,
								 ca_integrator,
								 ip3_integrator);
				     cell->compute_gap_junctions(dx_);
				}
			   }
		    )
	       );
     }
#else
     for(auto& c : cells_) {
	  c.do_step<integrator_t>(t_, dt, 
				  JPLCago,
				  ca_integrator,
				  ip3_integrator);
	  c.compute_gap_junctions(dx_);
     }
#endif // MULTI_THREAD
     post_update_step_(dt);
}


void halidi::System::print_system(output::Output_Streams& out) const
{
     const size_type SIZE(cells_.size());
     size_type points_printed(0);
     for (size_type i(0) ; i < SIZE ; ++i) {
	  out.v_out << t_.value() << " " << i << " ";
	  out.w_out << t_.value() << " " << i << " ";
	  cell_at_(i).print_cell(out, t_, points_printed);
     }
     out << std::endl;
}

void halidi::System::print_cell_at(size_type cell_idx,
				   size_type point_idx, 
				   output::Output_Streams& out) const
{
     out << t_.value()  << " ";
     cells_[cell_idx].print_point(out, point_idx);
     out << std::endl;
}

void halidi::System::print_cell_at(std::vector<Point> array,
				   size_type point_idx,
				   output::Output_Streams& out) const
{
     out << t_.value() << " ";
     for(const auto p  : array) {
	  cells_[p].print_point(out, point_idx);
     }
     out << std::endl;
}

void halidi::System::print_fluxes_at(size_type cell_idx,
				     size_type point_idx,
				     std::ostream& out) const
{
     cells_[cell_idx].print_fluxes(out, t_, point_idx);
     out << std::endl;
}

void halidi::System::print_fluxes_at(std::vector<Point> array,
				     size_type point_idx,
				     output::Flux_Streams& out) const
{
     for (const auto p : array) {
	  cells_[p].print_fluxes(out.get_stream(), t_, point_idx);
     }
     out << std::endl;
}

cellsim::concentration_ut halidi::System::get_c_at(size_type i) const
{
     return cell_at_(i).c_avg();
}

cellsim::concentration_ut halidi::System::get_I_at(size_type i) const
{
     return cell_at_(i).I_avg();
}

// ==============================================================================

halidi::Cell& halidi::System::cell_at_(size_type i)
{
#ifdef WITH_RANGE_CHECK
     range_check_(x, y);
#endif // WITH_RANGE_CHECK

     /*
      * we store the cells in the shape of a 2D matrix which is
      * represented by a single vector internally
      */
     return cells_[i];
}

halidi::Cell halidi::System::cell_at_(size_type i) const
{
#ifdef WITH_RANGE_CHECK
     range_check_(x, y);
#endif // WITH_RANGE_CHECK

     /*
      * we store the cells in the shape of a 2D matrix which is
      * represented by a single vector internally
      */
     return cells_[i];
}

void halidi::System::set_neighbours_(bool periodic_bc)
{     
     if (cells_.size() == 1) {
	  if (periodic_bc) {
	       cell_at_(0).set_neighbours(&cell_at_(0), &cell_at_(0));
	  }
     }
     else {
	  const auto last(end(cells_));
	  for (auto it(begin(cells_)), next(begin(cells_)+1), prev(last-1) ;
	       it != last ; ++it, ++next, ++prev) {
	       if (next == last) {
		    next = begin(cells_);
	       }
	       if (prev == last) {
		    prev = begin(cells_);
	       }
	       
	       it->set_neighbours(&(*prev), &(*next));
	  }
	  
	  // If not periodic BC, then unset neighbours where appropriate
	  if (!periodic_bc) {
	       cells_.front().set_neighbours(nullptr, &cell_at_(1));
	       auto it(end(cells_));
	       it -= 2;
	       cells_.back().set_neighbours(&(*it), nullptr);
	  }
     }
}

void halidi::System::post_update_step_(time_ut dt)
{
     /*
      * Once all cells have been taken care of, we can replace
      * the values of the current time with the ones from the next.
      */
     for(auto& c : cells_) {
	  c.swap_times();
     }

     t_ += dt;
}
