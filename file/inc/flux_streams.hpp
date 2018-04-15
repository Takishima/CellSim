/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef FLUX_STREAMS_HPP_INCLUDED
#define FLUX_STREAMS_HPP_INCLUDED

#include "definitions.hpp"

#include <fstream>
#include <string>
#include <vector>

#include <memory>

namespace std {
     template<typename T, typename ...Args>
     std::unique_ptr<T> make_unique( Args&& ...args )
     {
	  return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
     }
}

namespace cellsim {
namespace output {

     /*!
      *
      *
      */
     class Flux_Streams
     {
	  typedef std::basic_ostream<char, std::char_traits<char> > cout_t;
	  typedef cout_t& (*StandardEndLine)(cout_t&);
	  
	  typedef std::vector<std::unique_ptr<std::ofstream>> ofstream_vector_t;
	  // typedef std::vector<std::ofstream> ofstream_vector_t;
	  typedef ofstream_vector_t::iterator iterator_t;
	  
     public:
	  /*!
	   * \brief Simple constructor
	   *
	   * This class purpose is to simplify the creation and outputting of
	   * fluxes values for a selection of cells.
	   *
	   * \param ncells Number of cells to prepare files for
	   * \param output_file Name of the output file (used as a pattern here)
	   */	  
	  Flux_Streams(size_type ncells, std::string output_file);

	  /*!
	   * \brief Changes the precision of the output
	   *
	   * \param p New precision value
	   */
	  void precision(size_type p);

	  /*!
	   * \brief Output stream operator for standard IO manipulator
	   *
	   * \param manip Manipulator
	   * \return This object modified
	   * \note \c manip is passed to \e all the file streams.
	   */
	  Flux_Streams& operator<<(StandardEndLine manip);
	  
	  /*!
	   * \brief Returns an output stream
	   *
	   * Calling this method multiple time will cycle through all the 
	   * output streams stored in the Flux_Streams.
	   *
	   * \return Output stream
	   */
	  std::ostream& get_stream();

	  /*!
	   * \brief Generic output stream operator
	   *
	   * Outputs data to all streams
	   *
	   * \param t Value to output to all streams
	   * \return This object modified
	   */
	  template <typename T>
	  Flux_Streams& operator<<(T&& t)
	       {
		    for (auto& out : cell_flux_out_) {
			 *out << t;
		    }
		    return *this;
	       }

     private:
	  //! Array of files (one for each cell)
	  ofstream_vector_t cell_flux_out_;
	  iterator_t current_; //!< Iterator to the current file stream
     };
} // namespace output
} // namespace cellsim

#endif //FLUX_STREAMS_HPP_INCLUDED
