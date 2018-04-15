/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef SYSTEM_KOENIGSBERGER_IMPL_HPP_INCLUDED
#define SYSTEM_KOENIGSBERGER_IMPL_HPP_INCLUDED

#include "utility.hpp"

CLANG_DIAG_OFF(sign-conversion)

//#  include "impl/system_koenigsberger_thread_impl.hpp"

template <typename time_t, typename activator_t>
void cellsim::model_koenigsberger::System::activate(time_t dt, 
						    activator_t& activator)
{
     using integrators::eval_if;
     using equation_type::is_same_eq;
	  
     assert(dt.value() > 0.);

#ifdef MULTI_THREAD
     size_type start{0}, step{0}, change{0};
     std::tie(step, change) = utility::threads_step(cells_.size(),
     						    n_threads());
     std::vector<std::thread> threads;

     const auto t(t_);
     const auto Nx(Nx_);
     const auto BEGIN(begin(cells_));
     const auto END(end(cells_));
     for (auto it(begin(cells_)) ; it != END ; it += step, start += step) {
     	  if (start == change) {
     	       --step;
     	  }
	  
	  threads.push_back(
	       std::thread(
		    [=, &activator] () 
		    {
			 const auto last(it+step);
			 for(auto cell(it); cell != last ; ++cell) {
			      if (activator.do_activation((cell - BEGIN) % Nx)) {
				   cell->do_step<integrators::Forward_Euler>(t, dt, activator);
				   /*
				    * If activator is designed to act directly on a cell,
				    * apply it now that integration is done
				    */
				   eval_if<std::is_same<typename activator_t::equation_type,
							equation_type::NO_EQN>::value>::
					apply_activator(activator, *it, dt);

			      }
			      else {
				   cell->do_step<integrators::Forward_Euler>(t, dt);
			      }
			 }
		    }
		    )
	       );
     }
     
     // Wait for all threads to finish
     std::for_each(begin(threads), end(threads), [](std::thread& th) {th.join();});

#else
     const auto BEGIN(begin(cells_));
     const auto END(end(cells_));
     for(auto it(BEGIN) ; it != END ; ++it) {
	  if (activator.do_activation((it - BEGIN) % Nx_)) {
     	       it->do_step<integrators::Forward_Euler>(t_, dt, activator);

	       /*
		* If activator is designed to act directly on a cell,
		* apply it now that integration is done
		*/
	       eval_if<std::is_same<typename activator_t::equation_type,
				    equation_type::NO_EQN>::value>::
		    apply_activator(activator, *it, dt);
	  }
     	  else {
     	       it->do_step<integrators::Forward_Euler>(t_, dt);
     	  }
     }
#endif // MULTI_THREAD

     end_step_(dt);
}

CLANG_DIAG_ON(sign-conversion)

#endif //SYSTEM_KOENIGSBERGER_IMPL_HPP_INCLUDED
