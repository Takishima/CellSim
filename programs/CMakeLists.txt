include_directories( inc )
include_directories( ../core/inc )
include_directories( ../file/inc )
include_directories( ../integrators/inc )
include_directories( ../model/inc )

add_library( programs STATIC
  src/halidi.cpp
  src/koenigsberger.cpp
  src/halidi_bifurcation_diagram.cpp
  src/koenigs_bifurcation_diagram.cpp
 )