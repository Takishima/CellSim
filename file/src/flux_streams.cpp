#include "flux_streams.hpp"
#include "utility.hpp"


#include <iostream>


namespace output = cellsim::output;

output::Flux_Streams::Flux_Streams(size_type ncells, std::string output_file)
     : cell_flux_out_{}, current_{}
{
     const auto split = utility::split_filename(output_file);
     
     for (size_type i(0) ; i < ncells ; ++i) {
	  // cell_flux_out_.emplace_back(
	  //      utility::create_filename(split, 
	  // 				"flux", 
	  // 				utility::padded_str(i)
	  // 	    )
	  //      );
	  cell_flux_out_.push_back(
	       std::make_unique<std::ofstream>(
		    utility::create_filename(split, 
					     "flux", 
					     utility::padded_str(i)
			 )
		    )
	       );
     }
     current_ = begin(cell_flux_out_);
}

CLANG_DIAG_OFF(sign-conversion)

void output::Flux_Streams::precision(size_type p)
{
     for (auto& el : cell_flux_out_) {
	  el->precision(p);
	  // el.precision(p);
     }
}

CLANG_DIAG_ON(sign-conversion)

output::Flux_Streams& output::Flux_Streams::operator<<(StandardEndLine manip)
{
     for (auto& el : cell_flux_out_) {
	  // manip(el);
	  manip(*el);
     }
     return *this;
}

std::ostream& output::Flux_Streams::get_stream()
{
     auto& tmp = *current_;
     // auto tmp = current_;
     ++current_;
     if (current_ == end(cell_flux_out_)) {
	  current_ = begin(cell_flux_out_);
     }
     return *tmp;
}
