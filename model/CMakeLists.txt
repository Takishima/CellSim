include_directories( inc )
include_directories( ../core/inc )
include_directories( ../integrators/inc )

add_library( model STATIC
  src/constants_koenigsberger.cpp
  src/model_functions_2009.cpp
  src/model_equations_2009.cpp
  src/model_functions_2012.cpp
  src/model_equations_2012.cpp
  src/model_functions_koenigsberger.cpp
  src/model_equations_koenigsberger.cpp
 )