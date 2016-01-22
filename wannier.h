#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <armadillo>

#include "lattice.h"
#include "bloch.h"



#ifndef WANNIER_H_
#define WANNIER_H_


class Twannier{
public:

   	Twannier ( Tbloch bl1, int which_s );
    void setupXXmatrix();
    std::vector< std::complex<long double> > diagonalizeXXmatrix();
    void calc_wannier_f( );
    void rotate_to_real();
	void shift_to_real();
    void normalize_wannier();
    
    void print_wannier( std::string filename );

	double calc_energy();

//    void write_coeffs();
//	long double overlap_with( const Tmps_state& s1);
//    void normalize();

 //   std::vector< std::vector< std::complex<long double> > > return_bloch_functions(){ return bloch_functions; };
    
private:
	int which_site;
    int ksize;
    int lsize;
    double x_max;
    int n_of_bands;   

    Tlattice lattice;     
    
	std::vector< std::complex<long double> > eigVec;
    arma::Mat< std::complex<double> > XXm;

	std::vector< std::vector< double > > energies;
    std::vector< std::vector< std::vector< std::complex<long double> > > > bloch_functions;

    std::vector< std::pair < double,  std::complex<long double>  > > wannier_fun;
//	std::size_t length;
//	std::vector< std::complex<long double> > fourier_coeffs;
};

#endif


