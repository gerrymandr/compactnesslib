cmake_minimum_required (VERSION 3.9)

project (compactnesslib 
  VERSION 1.0.0
  DESCRIPTION "Measures compactness of shapes, especially political districts"
  LANGUAGES CXX
)

# add_subdirectory(apps)
# add_subdirectory(tests)

find_package(OpenMP)

add_library(compactnesslib
  src/bounded_scores.cpp
  src/csv.cpp
  src/geojson.cpp
  src/geom.cpp
  src/neighbours.cpp
  src/Props.cpp
  src/shapefile.cpp
  src/SpIndex.cpp
  src/unbounded_scores.cpp
  src/wkt.cpp

  lib/shapelib/dbfopen.cpp
  lib/shapelib/safileio.cpp
  lib/shapelib/shpopen.cpp
  lib/shapelib/shptree.cpp

  lib/clipper/clipper.cpp

  lib/spatialindex_c/rtree/BulkLoader.cpp
  lib/spatialindex_c/rtree/Index.cpp
  lib/spatialindex_c/rtree/Leaf.cpp
  lib/spatialindex_c/rtree/Node.cpp
  lib/spatialindex_c/rtree/RTree.cpp
  lib/spatialindex_c/rtree/Statistics.cpp
  lib/spatialindex_c/spatialindex/LineSegment.cpp
  lib/spatialindex_c/spatialindex/Point.cpp
  lib/spatialindex_c/spatialindex/Region.cpp
  lib/spatialindex_c/storagemanager/MemoryStorageManager.cpp
  lib/spatialindex_c/tools/Tools.cpp
)

target_include_directories(compactnesslib
PUBLIC
  api
  lib
)

target_compile_features(compactnesslib PUBLIC cxx_std_11)

target_link_libraries(compactnesslib
PRIVATE
  OpenMP::OpenMP_CXX
)

add_executable(unittest
  unittest/test.cpp
)

target_link_libraries(unittest
PRIVATE
  compactnesslib
)

target_compile_features(unittest PUBLIC cxx_std_11)

enable_testing()
add_test(NAME UnitTests COMMAND unittest WORKING_DIRECTORY unittest)
