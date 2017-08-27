SOURCES = $(wildcard ../*.cpp) $(wildcard ../shapelib/*.cpp) main.cpp Timer.cpp
OBJECTS = $(SOURCES:.cpp=.o)

CXX_FLAGS = --std=c++11 -march=native -mtune=native -g -O3 -DDOCTEST_CONFIG_DISABLE

all: $(OBJECTS)
	$(CXX) $(CXX_FLAGS) $(OBJECTS) -o compactness.exe -Wall -Wpedantic

%.o: %.cpp
	$(CXX) -c $(CXX_FLAGS) $< -o $@ 

clean:
	rm -rf compactness.exe ../*.o ../shapelib/*.o