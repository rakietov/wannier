#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>



#ifndef LATTICE_H_
#define LATTICE_H_

class Tlattice{
public:

    Tlattice(){ length = 0; };
   	Tlattice ( std::vector< std::complex<long double> > coeffs );
    ~Tlattice(){ };

    void write_coeffs();
    std::complex<long double> get_ith_fourier_coeff( int i ){ if(i < length) return fourier_coeffs[i]; else return (0.,0.); };
//	long double overlap_with( const Tmps_state& s1);
//    void normalize();

private:
	std::size_t length;
	std::vector< std::complex<long double> > fourier_coeffs;
};


#endif 
