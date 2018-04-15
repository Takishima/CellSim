/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef OUTPUT_STREAMS_HPP_INCLUDED
#define OUTPUT_STREAMS_HPP_INCLUDED

#include "definitions.hpp"

#include <fstream>
#include <string>

namespace cellsim {
namespace output {

     /*!
      * \brief Class to ease the output of all dynamical variables.
      *
      */
     class Output_Streams
     {
	  typedef std::basic_ostream<char, std::char_traits<char> > cout_t;
	  typedef cout_t& (*StandardEndLine)(cout_t&);

     public:
	  /*!
	   * \brief Simple constructor
	   *
	   * \param output_file Name of the output file (used as a pattern here)
	   * \param prefix String to prepend to output file name
	   */	  
	  Output_Streams(std::string output_file, std::string prefix = "");

	  /*!
	   * \brief Changes the precision of the output
	   *
	   * \param p New precision value
	   */
	  void precision(size_type p);

	  /*!
	   * \brief Changes the width of the output
	   *
	   * \param w New width value
	   */
	  void width(size_type w);

	  /*!
	   * \brief Output stream operator for standard IO manipulator
	   *
	   * \param manip Manipulator
	   * \return This object modified
	   */
	  Output_Streams& operator<<(StandardEndLine manip);

	  /*!
	   * \brief Output stream operator for values
	   *
	   * \param t Value to output
	   * \return This object modified
	   */
	  template <typename T>
	  Output_Streams& operator<<(T&& t)
	       {
		    ca_out << t;
		    ip3_out << t;
		    s_out << t;
		    v_out << t;
		    w_out << t;
		    return *this;
	       }

	  std::ofstream ca_out;
	  std::ofstream ip3_out;
	  std::ofstream s_out;
	  std::ofstream v_out;
	  std::ofstream w_out;
     };

} // namespace output
} // namespace cellsim

#endif //OUTPUT_STREAMS_HPP_INCLUDED
