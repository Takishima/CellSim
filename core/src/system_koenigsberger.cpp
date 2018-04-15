#include "system_koenigsberger.hpp"
#include "PLC_activator_koenigsberger.hpp"
#include "output_streams.hpp"
#include "flux_streams.hpp"

#ifdef MULTI_THREAD
#include <functional>
#endif // MULTI_THREAD

CLANG_DIAG_OFF(sign-conversion)

namespace model = cellsim::model_koenigsberger;

using namespace boost::units::si;

model::System::System(size_type Nx, size_type Ny, Constants_f_param_t constants)
     : System(Nx, Ny,
	      0._uM, 
	      0._uM,
	      0._uM,
	      0._mV,
	      0.,
	      constants)
{}

model::System::System(size_type Nx, 
		      size_type Ny,
		      concentration_ut c, 
		      concentration_ut I,
		      concentration_ut s, 
		      electric_potential_ut v, 
		      double w, 
		      Constants_f_param_t constants)
     : Nx_(Nx), Ny_(Ny), 
       t_(0. * second),
       cells_(Nx_ * Ny_, Cell(c, I, s, v, w, constants))
{
     set_neighbours_();

     // for(size_t x(0) ; x < Nx_ ; ++x) {
     // 	  for(size_t y(0) ; y < Ny_ ; ++y) {
     // 	       std::cout << "(" << x << ", " << y << ") : \t"
     // 			 << cell_at_(x, y).neighbours().size() << "\t";

     // 	       for(auto ptr : cell_at_(x, y).neighbours()) {
     // 		    if (ptr == nullptr) {
     // 			 std::cout << R"(PROBLEM)";
     // 		    }
     // 		    else {
     // 			 ptr->neighbours().size();
     // 		    }
     // 	       }
     // 	       std::cout << std::endl;
     // 	  }
     // }
}

// =============================================================================
// Public methods
void model::System::do_step(time_ut dt)
{
     assert(dt.value() > 0.);

#ifdef MULTI_THREAD     
     size_type start{0}, step{0}, change{0};
     std::tie(step, change) = utility::threads_step(cells_.size(),
						    n_threads());
     std::vector<std::thread> threads;

     const auto t(t_);
     const auto END(end(cells_));
     for (auto it(begin(cells_)) ; it != END ; it += step, start += step) {
	  if (start == change) {
	       --step;
	  }
	  
	  threads.push_back(
	       std::thread(
		    [=] () 
		    {
			 const auto last(it+step);
			 for(auto cell(it); cell != last ; ++cell) {
			      cell->do_step<integrators::Forward_Euler>(t, dt);
			 }
		    }
		    )
	       );
     }
     
     // Wait for all threads to finish
     std::for_each(begin(threads), end(threads), [](std::thread& th) {th.join();});


#else

     for (auto& c : cells_) {
	  // c.do_step<integrators::Standard_RK4>(t_, dt);
	  c.do_step<integrators::Forward_Euler>(t_, dt);
     }

#endif // MULTI_THREAD

     end_step_(dt);
}

void model::System::print_cell_at(size_type x, size_type y,
				  std::ofstream& out, TYPE value) const
{
     out << t_.value()  << " ";
     cell_at_(x, y).print(out, value);
     out << std::endl;
}

void model::System::print_cell_at(size_type x, size_type y,
				  output::Output_Streams& out) const
{
     out << t_.value()  << " ";
     cell_at_(x, y).print(out);
     // cell_at_(x, y).print(out.ca_out, cellsim::CA);
     // cell_at_(x, y).print(out.ip3_out, cellsim::IP3);
     // cell_at_(x, y).print(out.s_out, cellsim::S);
     // cell_at_(x, y).print(out.v_out, cellsim::V);
     // cell_at_(x, y).print(out.w_out, cellsim::W);
     out << std::endl;
}

void model::System::print_cell_at(std::vector<Point> array,
				  std::ofstream& out, TYPE value) const
{
     out << t_.value() << " ";
     for(const auto& p  : array) {
	  cell_at_(std::get<0>(p), std::get<1>(p)).print(out, value);
     }
     out << std::endl;
}
void model::System::print_cell_at(std::vector<Point> array,
				  output::Output_Streams& out) const
{
     out << t_.value() << " ";
     for(const auto& p  : array) {
	  cell_at_(std::get<0>(p), std::get<1>(p)).print(out.ca_out, 
							 cellsim::CA);
	  cell_at_(std::get<0>(p), std::get<1>(p)).print(out.ip3_out, 
							 cellsim::IP3);
	  cell_at_(std::get<0>(p), std::get<1>(p)).print(out.s_out, 
							 cellsim::S);
	  cell_at_(std::get<0>(p), std::get<1>(p)).print(out.v_out,
							 cellsim::V);
	  cell_at_(std::get<0>(p), std::get<1>(p)).print(out.w_out,
							 cellsim::W);
     }
     out << std::endl;
}

void model::System::print_fluxes_at(Point p, std::ofstream& out) const
{
     cell_at_(std::get<0>(p), std::get<1>(p)).print_fluxes(out, t_);
     out << std::endl;
}

void model::System::print_fluxes_at(std::vector<Point> array,
				     output::Flux_Streams& out) const
{
     for (const auto p : array) {
	  cell_at_(std::get<0>(p), std::get<1>(p)).print_fluxes(out.get_stream(), t_);
     }
     out << std::endl;
}
cellsim::concentration_ut model::System::get_c_at(size_type x, size_type y) const
{
     return cell_at_(x, y).c();
}

cellsim::concentration_ut model::System::get_I_at(size_type x, size_type y) const
{
     return cell_at_(x, y).I();
}


// =============================================================================
// Private methods

void model::System::print_fluxes_helper_(size_type x, size_type y,
					 std::ofstream& out) const
{
     using namespace model_functions;
     cell_at_(x, y).print_fluxes(out, t_);
}

void model::System::end_step_(time_ut dt)
{
     /*
      * Once we have stepped all the cells in time, we can swap the values
      * of the current time with the ones from the next.
      */
     for (auto& c : cells_) {
	  c.swap_times();
     }

     t_ += dt;
}


model::Cell& model::System::cell_at_(size_type x, size_type y)
{
#ifdef WITH_RANGE_CHECK
     range_check_(x, y);
#endif // WITH_RANGE_CHECK

     /*
      * we store the cells in the shape of a 2D matrix which is
      * represented by a single vector internally
      */
     return cells_[x + Nx_ * y];
}

model::Cell model::System::cell_at_(size_type x, size_type y) const
{
#ifdef WITH_RANGE_CHECK
     range_check_(x, y);
#endif // WITH_RANGE_CHECK

     /*
      * we store the cells in the shape of a 2D matrix which is
      * represented by a single vector internally
      */
     return cells_[x + Nx_ * y];
}

const model::Cell* model::System::address_at_(size_type x, size_type y) const
{
#ifdef WITH_RANGE_CHECK
     range_check_(x, y);
#endif // WITH_RANGE_CHECK

     return &cells_[x + Nx_ * y];
}

void model::System::set_neighbours_()
{
     const size_type Nx_1{Nx_ - 1}, Ny_1{Ny_ - 1};

     Cell::neighbours_array neighbours;

     // set neighbours
     for(size_t x(0) ; x < Nx_ ; ++x) {
     	  for(size_t y(0) ; y < Ny_ ; ++y) {

	       /*
		* Each cell has 6 neighbours :
		*   - 4 in the usual sense for a 2D lattice
		*   - 2 extra due to the fact that each other column 
		*     is shifted upwards:
		*
		*   *         #         *            
		*        #         #                 *    #    #    #    *
		*   *         x         *     ->     
		*        #         #                 *    #    x    #    *
		*   *         #         *            
		*        *         *                 *    *    #    *    *
		*/
 
	       // usual 4 neighbours
	       if (x != 0) {
		    neighbours.push_back(address_at_(x-1, y));
	       }
	       if (x != Nx_1) {
		    neighbours.push_back(address_at_(x+1, y));
	       }
	       if (y != 0) {
		    neighbours.push_back(address_at_(x, y-1));
	       }
	       if (y != Ny_1) {
		    neighbours.push_back(address_at_(x, y+1));
	       }

	       /*
		* We shift every other column upwards => 2 extra neighbours
		* x even :         x odd :
		*     #            # # #
		*   # x #          # x #
		*   # # #            #
		*/
	       if (x % 2 == 0) {
		    if (x != 0 && y != 0) {
			 neighbours.push_back(address_at_(x-1, y-1));
		    }
		    if (x != Nx_1 && y != 0) {
			 neighbours.push_back(address_at_(x+1, y-1));
		    }
	       }
	       else {
		    if (x != 0 && y != Ny_1) {
			 neighbours.push_back(address_at_(x-1, y+1));			 
		    }
		    if (x != Nx_1 && y != Ny_1) {
			 neighbours.push_back(address_at_(x+1, y+1));			 
		    }
	       }

	       cell_at_(x, y).set_neighbours(neighbours);
	       neighbours.clear();
     	  }

	  /*
	   * We need to swap the increments because the extra nearest neighbours
	   * change when x is even or odd.
	   *
	   * x even :         x odd :
	   *     #            # # #
	   *   # x #          # x #
	   *   # # #            #
	   */
	  // std::swap(shift1, shift1_old);
	  // std::swap(shift2, shift2_old);
     }
}

#ifdef WITH_RANGE_CHECK
void model::System::range_check_(size_type x, size_type y) const
{
     if (x >= Nx_) {
	  throw std::out_of_range("Out of range error: System::range_check (x)");
     }
     else if (y >= Ny_) {
	  throw std::out_of_range("Out of range error: System::range_check (y)");
     }
}
#endif // WITH_RANGE_CHECK

CLANG_DIAG_ON(sign-conversion)
