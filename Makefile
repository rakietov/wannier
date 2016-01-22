wannier: main.cpp Makefile
	 g++ *.cpp -o wannier -O3 -std=c++11 -march=native -I ~/armadillo-6.400.3/include -DARMA_DONT_USE_WRAPPER -DARMA_USE_ARPACK -fopenmp -llapack -larpack 
 
   
