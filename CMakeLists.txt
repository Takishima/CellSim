cmake_minimum_required(VERSION 2.6)
project (Cellsim CXX)

# Project defines (can be modified)
#add_definitions( -DWITH_RANGE_CHECK )
#add_definitions( -DMULTI_THREAD )
#add_definitions( -DHALIDI_2009 )
add_definitions( -DHALIDI_2012 )
#add_definitions( -DNO_UNITS )
add_definitions( -DSIMPLE_PROGRESS_BAR )

#add_definitions( -DPLC_ACTIVATION )
#add_definitions( -DELECTRIC_ACTIVATION )
add_definitions( -DPOTASSIUM_CHORIDE_ACTIVATION )

# ==============================================================================
#                                                                              #
# Below this point, everything should be fine whatever the configuration       #
#                                                                              #
# ==============================================================================

# Compiler flags
set( CMAKE_CXX_FLAGS "-std=c++11" )

set( CLANG_CXXFLAGS "-Weverything -Wno-c++98-compat -Wno-missing-prototypes" )
set( CLANG_CXXFLAGS "${CLANG_CXXFLAGS} -stdlib=libc++" )

set ( GCC_CXXFLAGS "-Wall -Wextra -ansi -pedantic -Wshadow -Wfloat-equal" )
set ( GCC_CXXFLAGS "${GCC_CXXFLAGS} -Wctor-dtor-privacy -Wnon-virtual-dtor" )
set ( GCC_CXXFLAGS "${GCC_CXXFLAGS} -Woverloaded-virtual -Wmissing-include-dirs" )
set ( GCC_CXXFLAGS "${GCC_CXXFLAGS} -Wconversion -Wwrite-strings -Wold-style-cast" )
set ( GCC_CXXFLAGS "${GCC_CXXFLAGS} -Wpointer-arith -Wcast-qual -Wcast-align" )
set ( GCC_CXXFLAGS "${GCC_CXXFLAGS} -Wswitch-enum -Wundef -Wredundant-decls" )
set ( GCC_CXXFLAGS "${GCC_CXXFLAGS} -Wstrict-null-sentinel -Wwrite-strings" )
set ( GCC_CXXFLAGS "${GCC_CXXFLAGS} -Wunreachable-code -Wno-unknown-pragmas" )
set ( GCC_CXXFLAGS "${GCC_CXXFLAGS} -Weffc++" )

# ==============================================================================
# Check which compiler we are dealing with and check its version
execute_process( COMMAND ${CMAKE_CXX_COMPILER} -v
  ERROR_VARIABLE _compiler_output )

# ======================================

set( IS_GCC 0 )
set( IS_CLANG 0 )

if ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
  set( IS_GCC 1 )
  set( CMAKE_CXX_FLAGS "${GCC_CXXFLAGS} ${CMAKE_CXX_FLAGS}" )
  string(REGEX REPLACE ".* ([0-9]\\.[0-9]\\.[0-9]) .*" "\\1"
    COMPILER_VERSION ${_compiler_output})

elseif ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" )
  set( IS_CLANG 1 )
  string(REGEX REPLACE ".* ([0-9]\\.[0-9]\\.[0-9]) .*" "\\1"
    COMPILER_VERSION ${_compiler_output})

  if ( ${COMPILER_VERSION} VERSION_GREATER 3.2 )
    set( CLANG_CXXFLAGS "${CLANG_CXXFLAGS} -Wno-documentation-deprecated-sync" )
    set( CLANG_CXXFLAGS "${CLANG_CXXFLAGS} -Wno-documentation" )
  endif()
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_CXXFLAGS}" )

else ()
  message( FATAL_ERROR "Unkown compiler -> does not know which flag to use..." )
endif()

message(STATUS "C++ compiler version: ${COMPILER_VERSION} [${CMAKE_CXX_COMPILER}]")

# ======================================
# Set Release/Debug compile flags

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -DNDEBUG -DBOOST_DISABLE_ASSERTS -O3")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None, Debug, Release (default)."
      FORCE)
    message(STATUS "Setting build type to: ${CMAKE_BUILD_TYPE}")
ENDIF(NOT CMAKE_BUILD_TYPE)

# ======================================

if( IS_GCC AND COMPILER_VERSION VERSION_LESS "4.7")
  message( FATAL_ERROR "This project requires at least GCC-4.7.x" )

elseif( IS_CLANG AND COMPILER_VERSION VERSION_LESS "3.1")
  message( FATAL_ERROR "This project requires at least clang-3.1 (got: ${COMPILER_VERSION})" )
endif()

# ==============================================================================
# Find out which git branch we are in
execute_process( COMMAND git rev-parse --abbrev-ref HEAD
  OUTPUT_VARIABLE git_branch_name 
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# ======================================
# Compute executable name
set( exec_name "cellsim-${git_branch_name}.x" )
message(STATUS "Executable name: ${exec_name}" )

add_executable( "${exec_name}" programs/src/main.cpp )

# ==============================================================================
# Add sub-directories and include paths
set( DIR_LIST activators core file integrators model programs )
foreach(d ${DIR_LIST} )
  include_directories(${Cellsim_SOURCE_DIR}/${d}/inc)
  add_subdirectory( ${Cellsim_SOURCE_DIR}/${d} )
endforeach(d)


# ======================================

#set(Boost_USE_STATIC_LIBS ON)
find_package( Boost 1.46.0 COMPONENTS program_options )
target_link_libraries("${exec_name}" ${Boost_LIBRARIES})
message( STATUS "Using system Boost::Program_Options: ${Boost_LIBRARIES}" )
include_directories(${Boost_INCLUDE_DIRS})
message( STATUS "Using Boost headers from ${Boost_INCLUDE_DIRS}" )

# ==============================================================================
# Add doxygen target
find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(DOC_OUTPUT_DIR "documentation")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in 
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Doxygen" VERBATIM
    )
  add_custom_target(clean_doc
    rm -rf ${DOC_OUTPUT_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Cleaning documentation" VERBATIM
    )
endif(DOXYGEN_FOUND)

# ==============================================================================

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)
message(STATUS "Executable output path: ${EXECUTABLE_OUTPUT_PATH}" )

# Beware of interdependencies
target_link_libraries( "${exec_name}" core file model programs activators core file model )

