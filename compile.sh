g++ -I. -c linear_int.cpp
g++ -std=c++11 -I. -fopenmp 2d-R2.cpp linear_int.o -lgsl -lgslcblas -o R2.exe
g++ -std=c++11 -I. -fopenmp 2d-R3.cpp linear_int.o -lgsl -lgslcblas -o R3.exe

