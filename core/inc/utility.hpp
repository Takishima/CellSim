/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef UTILITY_HPP_INCLUDED
#define UTILITY_HPP_INCLUDED

#include <string>
#include <tuple>

#include "definitions.hpp"

namespace cellsim {
namespace utility {
     /*!
      * \brief Extracts the path from a given string.
      *
      * \param str String containing a filename + path
      * \return Path to file
      * \note This function simply looks for the last '/' char in \c str
      * and returns everything before that position.
      */
     std::string extract_path(const std::string str);

     /*!
      * \brief Extracts the basename from a given string.
      *
      * \param str String containing a filename + path
      * \return Basename of the file
      * \note This function simply looks for the last '/' char in \c str
      * and returns everything after that position.
      */
     std::string basename(const std::string str);

     //! Helper type
     typedef std::tuple<std::string, std::string, std::string> split_file_t;

     /*!
      * \brief Splits a filename between path, basename and extension
      *
      * \param name Filename with path and extension
      * \return Tuple containing <path, basename, extension>
      *
      * \note The criteria to find the extension is simply to look for the last
      *       '.' of \c name if it exists.\n
      */
     split_file_t split_filename(const std::string name);

     /*!
      * \brief Creates a valid filename
      *
      * \param filename Splitted filename pattern
      * \param prefix String to prepend to the filename
      * \param suffix String to prepend to the filename
      * \return Filename
      * \sa split_filename
      */
     std::string create_filename(split_file_t filename, 
				 std::string prefix,
				 std::string suffix);


     /*!
      * \brief Padding numbers while converting to a string
      *
      * \param n Number to convert to string
      * \return String with a representation of \c n padded with zeroes
      * \note Currently padding for a minimum total length of 4
      */
     std::string padded_str(size_type n);

     /*!
      * \brief Convert EXEC_MODE to a string
      * 
      * \param mode EXEC_MODE value to convert
      * \return String representation of \c mode
      */
     std::string to_string(EXEC_MODE mode);

#ifdef MULTI_THREAD
    /*!
      * \brief Helper functions when using multi-threading.
      *
      * Given a data size, this function computes the number of steps to use
      * to split \c data_size into \c n chunks.
      *
      * \param data_size Size of data
      * \param n Final number of chunks
      * \return Tuple with <step, change>\n
      *         \c step is the index step to use initially
      *         \c change is the index at which the \c step value should be
      *                   decreased by 1.
      *
      * Example: \c data_size = 18 \c n = 4\n
      * Return value : <5, 10>\n
      * Data chunks: 0-5 5-10 10-14 14-18
      *    step =      5   5 -> 4     4
      */
     template <typename size_type>
     std::tuple<size_type, size_type>
     threads_step(const size_type data_size, const size_type n)
     {
	  if (data_size <= n) {
	       return std::make_tuple(1, data_size);
	  }
	  else {
	       size_type step(data_size / n);
	       const size_type remainder(data_size % n);
	       if (remainder != 0) {
		    ++step;
		    return std::make_tuple(step, step*remainder);
	       }
	       else {
		    return std::make_tuple(step, data_size);
	       }
	  }
     }
#endif // MULTI_THREAD

} // namespace utility
} // namespace cellsim

#endif //UTILITY_HPP_INCLUDED
