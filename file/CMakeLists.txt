include_directories( inc )
include_directories( ../core/inc )
include_directories( ../integrators/inc )
include_directories( ../model/inc )
include_directories( ../programs/inc )

add_library( file STATIC
  src/argument_parsing.cpp
  src/output_streams.cpp
  src/biffdiag_input_options.cpp
  src/read_constants.cpp
  src/flux_streams.cpp
  src/read_constants_halidi.cpp
  src/read_constants_koenigsberger.cpp
  src/generic_input_options.cpp
  src/read_parameters.cpp
  src/halidi_input_options.cpp
  src/simulation_input_options.cpp
  src/koenigsberger_input_options.cpp
 )