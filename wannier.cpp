#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>


#include <armadillo>
#include "wannier.h"
#include "bloch.h"
#include "lattice.h"



// ----------------------------------------------------------------------------------------
Twannier::Twannier( Tbloch bl1, int which_s)
{
	which_site = which_s;
    ksize = bl1.return_ksize();
    lsize = bl1.return_lsize();
    x_max = bl1.return_x_max();
    n_of_bands = bl1.return_n_of_bands();  
    lattice = bl1.return_lattice();    
    bloch_functions = bl1.return_bloch_functions();

	this -> setupXXmatrix();
	eigVec =  diagonalizeXXmatrix();

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Twannier::setupXXmatrix()
{

    //std::cout << "setup XX matrix" << std::endl;
    XXm.resize( (2*ksize + 1)*n_of_bands, (2*ksize + 1)*n_of_bands);


	for(int i0 = 0; i0 < n_of_bands; ++i0)
	{
		for(int i0p = 0; i0p < n_of_bands; ++i0p)
		{
			for(int i1 = - ksize; i1 < ksize+1; ++i1)
			{		
				for( int i2 = -ksize; i2 < ksize+1; ++i2)
				{
					std::complex<long double> mat_ele = (0.,0.);
					for(int n1 = 0; n1 <= 2*lsize; ++n1)
					{
						for(int n2 = 0; n2 <= 2*lsize ; ++n2)
						{
							if( i1 - i2 + (2*ksize + 1)*(n1-n2) != 0)
							{
								std::complex<long double> tmp =  std::complex<long double>(
										double(std::pow((-1.), i1-i2+n1-n2)) / 
										((double(i1-i2)/(2.*ksize+1.)+ n1-n2 )));
								tmp *= std::complex<long double>(0.,0.5)  *
									   std::conj( bloch_functions[i0p][i2+ksize][n2] ) *  
										bloch_functions[i0][i1+ksize][n1];
								mat_ele += tmp;
							}
						}
					}	
					XXm(i1+ksize + i0*(2*ksize + 1),i2+ksize + i0p*(2*ksize+1) ) = mat_ele;
		            //std::cout << i1+ksize + i0p*(2*ksize + 1)<<" "<<i2+ksize + i0*(2*ksize+1)<<" "<<
					//		XXm(i1+ksize + i0p*(2*ksize + 1),i2+ksize + i0*(2*ksize+1) ) << std::endl;
				}
			}
		}
	}
    //std::cout << "setup XX matrix DONE" << std::endl;
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
std::vector< std::complex<long double> > Twannier::diagonalizeXXmatrix()
{
	arma::Mat< std::complex<double> > eigvec;
	arma::vec eigval;

	arma::eig_sym( eigval, eigvec, XXm);

	for(int i = 0; i < (2*ksize+1)*n_of_bands; ++i)
	{
		std::cout.precision(15);
		std::cout << eigval( i ) << std::endl;
	}
	for(int i = 0; i < (2*ksize+1)*n_of_bands-1; ++i)
	{
		std::cout.precision(15);
		std::cout << eigval( i+1 ) - eigval(i) << std::endl;
	}

	std::vector< std::complex<long double> > eigV;
    for(int i2 = 0; i2 < (2*ksize+1)*n_of_bands; i2++)
    {
        eigV.push_back( eigvec( i2, (2*ksize+1) + which_site )  );
    }
   
//	std::cout << __func__;  

	return eigV;

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Twannier::print_wannier( std::string fname )
{

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Twannier::shift_to_real()
{

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Twannier::rotate_to_real()
{

}
// ========================================================================================



// ----------------------------------------------------------------------------------------
void Twannier::normalize_wannier()
{

}
// ========================================================================================



// ----------------------------------------------------------------------------------------
void Twannier::calc_wannier_f()
{


    for(int ix = 0; ix < 10000; ix++)
    {
        double x = ix/100. - 50.;
        std::complex<long double> wannier_x (0., 0.);
		for(int i0 = 0; i0 < n_of_bands; ++i0)
		{
	        for(int i1= -ksize; i1 < ksize +1; ++i1)
		    {
				std::complex<long double> tmp = 0.;
				for(int i2 = -lsize; i2 < lsize+1; ++i2)
				{
	                tmp += bloch_functions[i0][i1+ksize][i2+lsize] * std::exp( std::complex<long double>
							(0., (i2 + i1 *1./(2.*ksize+1.))* x)   );
	                //std::cout << tmp << std::endl;
		        }
		        wannier_x += eigVec[i1+ksize + (2*ksize+1)*i0] * tmp;
			}
			std::cout << x << " " << wannier_x.real()<<" "<<wannier_x.imag()<< std::endl;
		}
	}
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
double Twannier::calc_energy()
{

}
// ========================================================================================









