#include "utility.hpp"

CLANG_DIAG_OFF(sign-conversion)

namespace utility = cellsim::utility;

cellsim::size_type helper(std::string str)
{
     return str.find_last_of("/");
}

std::string utility::extract_path(const std::string str)
{
     const size_type idx(helper(str));
     if (idx == std::string::npos) {
	  return "./";
     }
     else {
	  return std::string(str.begin(), str.begin()+idx);
     }
}

std::string utility::basename(const std::string str)
{
     const size_type idx(helper(str));
     if (idx == std::string::npos) {
	  return str;
     }
     else {
	  return std::string(str.begin()+idx, str.end());
     }     
}

std::tuple<std::string, std::string, std::string>
utility::split_filename(std::string name)
{
     std::string path;
     size_t idx(helper(name));
     if (idx == std::string::npos) {
	  path = "./";
     }
     else {
	  ++idx;
	  path = std::string(begin(name), begin(name)+idx);
	  name = name.substr(idx);
     }

     idx = name.find_last_of(".");

     if (idx == std::string::npos) {
	  return std::make_tuple(path, name, "");
     }
     else {
	  return std::make_tuple(
	       path,
	       std::string(name.begin(), name.begin() + idx),
	       std::string(name.begin() + idx, name.end())
	       );
     }
}

std::string utility::create_filename(split_file_t filename, 
				     std::string prefix,
				     std::string suffix)
{
     if (!std::get<1>(filename).empty()) {
	  if (!prefix.empty()) {
	       prefix += "_";
	  }
	  if (!suffix.empty()) {
	       suffix = "_" + suffix;
	  }
     }
     else if (!prefix.empty()) {
	  prefix += "_";
     }
     return std::get<0>(filename) 
	  + prefix + std::get<1>(filename) + suffix
	  + std::get<2>(filename);
}

std::string utility::padded_str(size_type n)
{
     std::string result;
     if (n < 10) {
	  result += '0';
     }
     if (n < 100) {
	  result += '0';
     }
     if (n < 1000) {
	  result += '0';
     }
     return (result += std::to_string(n));
}

std::string utility::to_string(EXEC_MODE mode)
{
     if (mode == BIFFDIAG_HALIDI) {
	  return "BIFFDIAG_HALIDI";
     }
     else if (mode == HALIDI) {
	  return "HALIDI";
     }
     else if (mode == BIFFDIAG_KOENIGS) {
	  return "BIFFDIAG_KOENIGS";
     }
     else if (mode == KOENIGSBERGER) {
	  return "KOENIGSBERGER";
     }
     else if (mode == INVALID) {
	  return "INVALID";
     }
     else {
	  return "unknown EXEC_MODE value";
     }
}

CLANG_DIAG_ON(sign-conversion)
