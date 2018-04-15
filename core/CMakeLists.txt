include_directories( inc )
include_directories( ../activators/inc )
include_directories( ../integrators/inc )
include_directories( ../file/inc )
include_directories( ../model/inc )

add_library( core STATIC
  src/point.cpp
  src/system_koenigsberger.cpp
  src/point_container.cpp
  src/utility.cpp
  src/cell.cpp
  src/point_index.cpp
  src/cell_koenigsberger.cpp
  src/system.cpp
  src/definitions.cpp
 )