#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>

//#include </home/sierant/Eigen324/Eigen/Dense>
#include <armadillo>
#include "wannier.h"
#include "bloch.h"
#include "lattice.h"



// ----------------------------------------------------------------------------------------
Twannier::Twannier( Tbloch bl1, int which_s) : 
	which_site (which_s), ksize(bl1.return_ksize()), lsize(bl1.return_lsize() ), x_max(bl1.return_x_max()),
	n_of_bands( bl1.return_n_of_bands()), lattice( bl1.return_lattice()), 
	bloch_functions(bl1.return_bloch_functions()), energies( bl1.return_energies() )
{

	this -> setupXXmatrix();
	eigVec =  diagonalizeXXmatrix();
	this -> calc_wannier_f();
	this -> normalize_wannier();
	this -> shift_to_real();
	std::cout <<"Energy = " << calc_energy() << std::endl;
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
					for(int n1 = -lsize; n1 < lsize + 1; ++n1)
					{
						for(int n2 = -lsize; n2 < lsize + 1; ++n2)
						{
							if( i1 - i2 + (2*ksize + 1)*(n1-n2) != 0)
							{
								std::complex<long double> tmp =  std::complex<long double>(
										double(std::pow((-1.), i1-i2+n1-n2)) / 
										((double(i1-i2)/(2.*ksize+1.)+ double(n1-n2) )));
								tmp *= std::complex<long double>(0.,0.5)  *
									   std::conj( bloch_functions[i0p][i2+ksize][n2+lsize] ) *  
										bloch_functions[i0][i1+ksize][n1+lsize];
								mat_ele += tmp;
								//std::cout << tmp << " " << mat_ele  << std::endl;
							}
						}
					}	
					XXm(i1+ksize + i0*(2*ksize + 1),i2+ksize + i0p*(2*ksize+1) ) = mat_ele;
					//std::cout<< "mat ele " << mat_ele <<std::endl;
		            //std::cout << i1+ksize + i0*(2*ksize + 1)<<" "<<i2+ksize + i0p*(2*ksize+1)<<" "<<
					//		XXm(i1+ksize + i0*(2*ksize + 1),i2+ksize + i0p*(2*ksize+1) ) << std::endl;
				}
			}
		}
	}
    std::cout << "setup XX matrix DONE" << std::endl;

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
std::vector< std::complex<long double> > Twannier::diagonalizeXXmatrix()
{
	arma::Mat< std::complex<double> > eigvec;
	arma::vec eigval;
	arma::eig_sym( eigval, eigvec, XXm);
	

	std::cout<< __func__ <<std::endl;

	
	//    Eigen::SelfAdjointEigenSolver< Eigen::Matrix<std::complex< long double>, 
	//						  Eigen::Dynamic, Eigen::Dynamic> > eigensolver(XXm);
//	std::cout << eigensolver.eigenvalues() <<std::endl;


for(int i = 0; i < (2*ksize+1)*n_of_bands; ++i)
	{
		std::cout.precision(15);
		std::cout << eigval( i ) << std::endl;
	}
	for(int i = 0; i < (2*ksize+1)*n_of_bands-1; ++i)
	{
		std::cout.precision(15);
		std::cout <<eigval(i+1) - eigval(i) <<std::endl;
	}

	std::vector< std::complex<long double> > eigV;
    for(int i2 = 0; i2 < (2*ksize+1)*n_of_bands; ++i2)
    {
		eigV.push_back( eigvec( i2, n_of_bands*ksize + which_site) );   
		//eigV.push_back( eigensolver.eigenvectors().col( n_of_bands*ksize  + which_site )[i2]  );
    }
   
//	std::cout << __func__;  

	return eigV;

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Twannier::print_wannier( std::string fname )
{
	
	std::fstream fs;
	fs.open( fname.c_str(), std::fstream::out);
	fs.precision(15);
	for(int i = 0; i < wannier_fun.size(); ++i)
	{
		fs << wannier_fun[i].first << " " << (wannier_fun[i].second).real() << 
			" " << (wannier_fun[i].second).imag() << std::endl;

	}
	fs.close();

}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Twannier::normalize_wannier()
{
	long double norm_sq = 0.0;
	for(int i = 0; i < wannier_fun.size()-1; ++i)
	{
		norm_sq += wannier_fun[i].second.real() * wannier_fun[i].second.real() * 
				   (wannier_fun[i+1].first - wannier_fun[i].first);
	
		norm_sq += wannier_fun[i].second.imag() * wannier_fun[i].second.imag() * 
				   (wannier_fun[i+1].first - wannier_fun[i].first);

	}

	for(int i = 0; i < wannier_fun.size(); ++i)
	{
		wannier_fun[i].second /= std::sqrt(norm_sq);
	}


}
// ========================================================================================


// ----------------------------------------------------------------------------------------
void Twannier::rotate_to_real()
{

}
// ========================================================================================



// ----------------------------------------------------------------------------------------
void Twannier::shift_to_real()
{

	double norm_re = 0.;
	double norm_im = 0.;
	for(int i = 0; i < wannier_fun.size()-1; ++i)
	{
		norm_re += wannier_fun[i].second.real() * wannier_fun[i].second.real();
		norm_im += wannier_fun[i].second.imag() * wannier_fun[i].second.imag();
	}

	if(norm_re > norm_im)
	{	
		for(int i = 0; i < wannier_fun.size(); ++i)
		{
			wannier_fun[i].second = std::complex<long double>(
				std::sqrt(wannier_fun[i].second.real() * wannier_fun[i].second.real() +
				wannier_fun[i].second.imag() * wannier_fun[i].second.imag())*
			    wannier_fun[i].second.real()/std::fabs( wannier_fun[i].second.real() ),0.);
		}
	}
	else
	{
		for(int i = 0; i < wannier_fun.size(); ++i)
		{
			wannier_fun[i].second = std::complex<long double>(
				std::sqrt(wannier_fun[i].second.real() * wannier_fun[i].second.real() +
				wannier_fun[i].second.imag() * wannier_fun[i].second.imag())*
			    wannier_fun[i].second.imag()/std::fabs( wannier_fun[i].second.imag() ),0.);
		}
	}

}
// ========================================================================================



// ----------------------------------------------------------------------------------------
void Twannier::calc_wannier_f()
{
    for(int ix = 0; ix < 10000; ix++)
    {
        long double x = ix/100. - 50.;
        std::complex<long double> wannier_x (0., 0.);
		for(int i0 = 0; i0 < n_of_bands; ++i0)
		{
	        for(int i1= -ksize; i1 < ksize +1; ++i1)
		    {
				std::complex<long double> tmp = std::complex<long double>(0., 0.);
				for(int i2 = -lsize; i2 < lsize+1; ++i2)
				{
	                tmp += bloch_functions[i0][i1+ksize][i2+lsize] * std::exp( std::complex<long double>
							(0., (2.*i2 + i1 *2./(2.*ksize+1.))* x) );
	                //std::cout << bloch_functions[i0][i1+ksize][i2+lsize] *
					//				std::exp( std::complex<long double>      
                    //         (0., (2.*i2 + i1 *2./(2.*ksize+1.))* x) ) << std::endl;
		        }
		        wannier_x += eigVec[i1+ksize + (2*ksize+1)*i0] * tmp;
				//std::cout <<"d_wan " << eigVec[i1+ksize + (2*ksize+1)*i0] * tmp <<std::endl;
			}
		}
		wannier_fun.push_back( std::make_pair(x, wannier_x) );
	}
	//std::cout<< __func__<< " " << wannier_fun.size() << std::endl;
}
// ========================================================================================


// ----------------------------------------------------------------------------------------
double Twannier::calc_energy()
{
	double energ = 0.0;

	for(int i0 = 0; i0 < n_of_bands; ++i0)
	{
		for(int i1 = -ksize; i1 < ksize + 1; ++i1)
		{
			energ +=std::pow(std::abs( eigVec[i0*(2*ksize+1) + i1+ksize]  ),2.) * energies[i0][i1+ksize];
		}

	}
	return energ;


}
// ========================================================================================









