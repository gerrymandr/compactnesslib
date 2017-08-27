SOURCES = $(wildcard ../*.cpp) $(wildcard ../shapelib/*.cpp) test.cpp
OBJECTS = $(SOURCES:.cpp=.o)

CXX_FLAGS = --std=c++11 -march=native -mtune=native -g --coverage -O0

all: $(OBJECTS)
	$(CXX) $(CXX_FLAGS) $(OBJECTS) -o compactness_test.exe  -Wall -Wpedantic

%.o: ./%.cpp
	$(CXX) -c $(CXX_FLAGS) $< -o $@ 

clean:
	rm -rf compactness_test.exe ../*.o ../shapelib/*.o *.gcno *.gcda *.gcov