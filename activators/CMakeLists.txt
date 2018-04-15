include_directories( inc )
include_directories( ../core/inc )
include_directories( ../integrators/inc )
include_directories( ../model/inc )

add_library( activators STATIC
  src/halidi_activator_base.cpp
  src/koenigsberger_activator_base.cpp
  src/PLC_activator.cpp
  src/PLC_activator_koenigsberger.cpp
  src/potassium_chloride_activator.cpp
  src/potassium_chloride_activator_koenigsberger.cpp
  src/potential_activator_base.cpp
  src/potential_activator.cpp
  src/potential_activator_tanh.cpp
  src/potential_activator_koenigsberger.cpp
  src/potential_activator_koenigsberger_tanh.cpp
 )