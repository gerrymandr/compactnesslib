all:
	g++ -O3 -g -DDOCTEST_CONFIG_DISABLE --std=c++11 -fopenmp -march=native -mtune=native main.cpp Timer.cpp ../*cpp ../shapelib/*cpp -o compactness.exe -Wall -Wpedantic

test:
	g++ --coverage -O0 -g --std=c++11 ../*cpp ../shapelib/*cpp test.cpp -o compactness_test.exe  -Wall -Wpedantic

cleantest:
	rm -rf compactness_test.exe *.gcno *.gcda *.gcov