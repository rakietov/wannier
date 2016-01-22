#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include "lattice.h"


// ----------------------------------------------------------------------------------------
Tlattice::Tlattice( std::vector< std::complex<long double> > coeffs  ) : length( coeffs.size() ) 
{
	for( int i = 0; i < coeffs.size(); ++i)
	{
		fourier_coeffs.push_back( coeffs[i]) ;
	}
}
// ========================================================================================

// ----------------------------------------------------------------------------------------
void Tlattice::write_coeffs()
{
	for( int i = 0; i < fourier_coeffs.size(); ++i)
	{
		std::cout<<fourier_coeffs[i]<<std::endl;
	}
}
// ========================================================================================

