#include "output_streams.hpp"
#include "utility.hpp"

namespace output = cellsim::output;

output::Output_Streams::Output_Streams(std::string output_file, 
				       std::string prefix)
     : ca_out{}, ip3_out{}, s_out{}, v_out{}, w_out{}
{
     const auto split = utility::split_filename(output_file);

     ca_out.open(utility::create_filename(split, prefix, "ca"));
     ip3_out.open(utility::create_filename(split, prefix, "ip3"));
     s_out.open(utility::create_filename(split, prefix, "s"));
     v_out.open(utility::create_filename(split, prefix, "v"));
     w_out.open(utility::create_filename(split, prefix, "w"));
}

output::Output_Streams& output::Output_Streams::operator<<(StandardEndLine manip)
{
     manip(ca_out);
     manip(ip3_out);
     manip(s_out);
     manip(v_out);
     manip(w_out);
     return *this;
}

CLANG_DIAG_OFF(sign-conversion)

void output::Output_Streams::precision(size_type p)
{
     ca_out.precision(p);
     ip3_out.precision(p);
     s_out.precision(p);
     v_out.precision(p);
     w_out.precision(p);
}

void output::Output_Streams::width(size_type w)
{
     ca_out.width(w);
     ip3_out.width(w);
     s_out.width(w);
     v_out.width(w);
     w_out.width(w);
}

CLANG_DIAG_ON(sign-conversion)
