#include "definitions.hpp"

cellsim::size_type& cellsim::n_threads()
{
     static size_type n_threads_{4};
     return n_threads_;
}
