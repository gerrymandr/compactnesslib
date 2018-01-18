SOURCES = $(wildcard ../*.cpp) $(wildcard ../shapelib/*.cpp) $(wildcard ../lib/*.cpp) test.cpp ../lib/spatialindex_c/spatialindex/Region.cpp ../lib/spatialindex_c/rtree/RTree.cpp ../lib/spatialindex_c/spatialindex/Point.cpp ../lib/spatialindex_c/tools/Tools.cpp ../lib/spatialindex_c/storagemanager/MemoryStorageManager.cpp ../lib/spatialindex_c/spatialindex/LineSegment.cpp ../lib/spatialindex_c/rtree/Leaf.cpp ../lib/spatialindex_c/rtree/Node.cpp ../lib/spatialindex_c/rtree/Statistics.cpp ../lib/spatialindex_c/rtree/BulkLoader.cpp ../lib/spatialindex_c/rtree/Index.cpp
OBJECTS = $(SOURCES:.cpp=.o)

CXX_FLAGS = --std=c++11 -march=native -mtune=native -g --coverage -O0

all: $(OBJECTS)
	$(CXX) $(CXX_FLAGS) $(OBJECTS) -o compactness_test.exe  -Wall -Wpedantic

%.o: ./%.cpp
	$(CXX) -c $(CXX_FLAGS) $< -o $@ 

clean:
	rm -rf compactness_test.exe *.o ../*.o ../shapelib/*.o ../lib/*.o *.gcno *.gcda *.gcov