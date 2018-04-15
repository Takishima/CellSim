/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef SYSTEM_IMPL_HPP_INCLUDED
#define SYSTEM_IMPL_HPP_INCLUDED

#ifdef MULTI_THREAD
#  include <thread>
#endif // MULTI_THREAD

// // template <typename activator_t>
// void System::do_step(time_ut dt,
// 		     flux_ut JPLCago, 
// 		     const ,activator_t activator)
// {
//     using namespace cellsim::equation_type;
    
//     integrators::FTCS<CA_EQN> ca_integrator(dt, dx_);
//     integrators::FTCS<IP3_EQN> ip3_integrator(dt, dx_);

//     const size_type SIZE(cells_.size());

// #ifdef MULTI_THREAD
//     size_type start{0}, step{0}, change{0};
//     std::tie(step, change) = utility::threads_step(cells_.size(),
//                                                    System::n_threads);
//     std::vector<std::thread> threads;

// #else
//     for(size_type i(0) ; i < SIZE ; ++i) {
//         if (activator.do_activation(i)) {
//             cells_[i].do_step<integrators::Forward_Euler>(t_, dt, 
//                                                           JPLCago,
//                                                           ca_integrator,
//                                                           ip3_integrator,
//                                                           activator);
//         }
//         else {
//             cells_[i].do_step<integrators::Forward_Euler>(t_, dt, 
//                                                           JPLCago,
//                                                           ca_integrator,
//                                                           ip3_integrator);
//         }
//         cells_[i].compute_gap_junctions(dx_);

//         cells_[i].swap_times();      
//     }
    
// #endif // MULTI_THREAD
    
// }

#endif //SYSTEM_IMPL_HPP_INCLUDED
