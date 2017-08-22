#!/bin/bash
g++ -g -DDOCTEST_CONFIG_DISABLE --std=c++11 *.cpp ../*cpp ../shapelib/*c -fopenmp
