cmake_minimum_required (VERSION 3.9)

project (compactnesslib
  VERSION 1.0.0
  DESCRIPTION "Measures compactness of shapes, especially political districts"
  LANGUAGES CXX
)

# add_subdirectory(apps)
# add_subdirectory(tests)

find_package(OpenMP)

add_library(libspatialindex
  submodules/libspatialindex/src/mvrtree/Index.cc
  submodules/libspatialindex/src/mvrtree/Node.cc
  submodules/libspatialindex/src/mvrtree/Leaf.cc
  submodules/libspatialindex/src/mvrtree/Statistics.cc
  submodules/libspatialindex/src/mvrtree/MVRTree.cc
  submodules/libspatialindex/src/tprtree/Index.cc
  submodules/libspatialindex/src/tprtree/Node.cc
  submodules/libspatialindex/src/tprtree/Leaf.cc
  submodules/libspatialindex/src/tprtree/Statistics.cc
  submodules/libspatialindex/src/tprtree/TPRTree.cc
  submodules/libspatialindex/src/storagemanager/Buffer.cc
  submodules/libspatialindex/src/storagemanager/RandomEvictionsBuffer.cc
  submodules/libspatialindex/src/storagemanager/DiskStorageManager.cc
  submodules/libspatialindex/src/storagemanager/MemoryStorageManager.cc
  submodules/libspatialindex/src/rtree/RTree.cc
  submodules/libspatialindex/src/rtree/BulkLoader.cc
  submodules/libspatialindex/src/rtree/Index.cc
  submodules/libspatialindex/src/rtree/Node.cc
  submodules/libspatialindex/src/rtree/Leaf.cc
  submodules/libspatialindex/src/rtree/Statistics.cc
  submodules/libspatialindex/src/capi/IdVisitor.cc
  submodules/libspatialindex/src/capi/Utility.cc
  submodules/libspatialindex/src/capi/ObjVisitor.cc
  submodules/libspatialindex/src/capi/Index.cc
  submodules/libspatialindex/src/capi/DataStream.cc
  submodules/libspatialindex/src/capi/BoundsQuery.cc
  submodules/libspatialindex/src/capi/CustomStorage.cc
  submodules/libspatialindex/src/capi/sidx_api.cc
  submodules/libspatialindex/src/capi/CountVisitor.cc
  submodules/libspatialindex/src/capi/LeafQuery.cc
  submodules/libspatialindex/src/capi/Error.cc
  submodules/libspatialindex/src/tools/rand48.cc
  submodules/libspatialindex/src/tools/Tools.cc
  submodules/libspatialindex/src/spatialindex/Point.cc
  submodules/libspatialindex/src/spatialindex/LineSegment.cc
  submodules/libspatialindex/src/spatialindex/TimePoint.cc
  submodules/libspatialindex/src/spatialindex/MovingPoint.cc
  submodules/libspatialindex/src/spatialindex/TimeRegion.cc
  submodules/libspatialindex/src/spatialindex/SpatialIndexImpl.cc
  submodules/libspatialindex/src/spatialindex/MovingRegion.cc
  submodules/libspatialindex/src/spatialindex/Region.cc
)
target_include_directories(libspatialindex PUBLIC submodules/libspatialindex/include)


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
PUBLIC
  libspatialindex
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
