#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>

#include "lattice.h"

#ifndef BLOCH_H_
#define BLOCH_H_

class Tbloch{
public:

//    Tbloch(){ length = 0; };
   	Tbloch ( Tlattice latt, int ksize, int lsize, double x_max, int n_of_bands, int which_b );
//    ~Tbloch(){ };

    void write_coeffs();
	void print_dispersion( std::string fname );
//	long double overlap_with( const Tmps_state& s1);
//    void normalize();

    int return_ksize(){ return ksize; };
    int return_lsize(){ return lsize; };
    double return_x_max(){ return x_max; };
    int return_n_of_bands() { return n_of_bands; };
    Tlattice return_lattice(){ return lattice; };

    std::vector< std::vector< std::vector< std::complex<long double> > > > return_bloch_functions(){ 
		return bloch_functions; };
	std::vector<std::vector<double > > return_energies(){ return energies; };
    
private:
	int which_band;
    int ksize;
    int lsize;
    double x_max;
    int n_of_bands;   
    Tlattice lattice;     
    
	std::vector< std::vector< double > > energies;
    std::vector< std::vector< std::vector< std::complex<long double> > > > bloch_functions;
//	std::size_t length;
//	std::vector< std::complex<long double> > fourier_coeffs;
};


#endif 


