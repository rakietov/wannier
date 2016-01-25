#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
//#include <armadillo>

#include </home/sierant/Eigen324/Eigen/Dense>

#include "bloch.h"

// ----------------------------------------------------------------------------------------
Tbloch::Tbloch( Tlattice latt, int ks, int ls,  double x_m, int bands, int which_b  ) 
{
	which_band =  bands * which_b;
    ksize = ks;
    lsize = ls;
    x_max = x_m;
    n_of_bands = bands;  
    lattice = latt;

// setup the Hamiltonian 
	Eigen::Matrix<  std::complex<long double>, Eigen::Dynamic, Eigen::Dynamic> ham;
	//arma::Mat< std::complex<double> > ham;
	ham.resize( 2*lsize + 1, 2*lsize + 1);

    for(int i0 = 0; i0 < n_of_bands; ++i0)
    {
        std::vector< std::vector< std::complex<long double> > > tmp;
        bloch_functions.push_back( tmp );              
		std::vector< double  > tmp2;
		energies.push_back( tmp2 );
    }

    for(int i = - ksize; i < ksize+1; ++i)
    {
        double kbier = i * 2./(2.*ksize+1.);
           
        for(int i1 = -lsize; i1 < lsize+1; ++i1)
        {
            for(int i2 = -lsize; i2 < lsize+1; ++i2)
            {
                ham(i1 + lsize, i2 + lsize) = std::complex<double>(0., 0.);        
            }
        }

        // kinetic energy
        for(int i1 = -lsize; i1 <= lsize; ++i1)
        {
            ham(i1 + lsize, i1 + lsize) = (2.*i1+kbier)*(2.*i1+kbier);
        }        

        for(int i1 = -lsize; i1 <= lsize; ++i1)
        {
            for(int i2 = -lsize; i2 <= lsize; ++i2)
            {
                std::complex<long double> dh = lattice.get_ith_fourier_coeff( abs( i1 - i2)  );
                if(i1 < i2)dh = std::conj(dh);
                //std::cout << i1+lsize << " " << i2+lsize << " "<< dh << std::endl;
                ham(i1 + lsize, i2 + lsize) += dh;   
            }
        }

	    Eigen::SelfAdjointEigenSolver< Eigen::Matrix<std::complex<long double>, 
							  Eigen::Dynamic, Eigen::Dynamic> > eigensolver(ham);
	
		//arma::vec eigval;
		//arma::Mat< std::complex<double> > eigvec;

		//arma::eig_sym(eigval, eigvec, ham, "std");

        
		//std::cout << kbier << std::endl;
        //std::cout << eigensolver.eigenvalues() <<std::endl<<std::endl;  

        for(int i0 = 0; i0 < n_of_bands; ++i0)
        {
            std::vector< std::complex<long double> > eigV;        
            for(int i2 = 0; i2 < 2*lsize+1; i2++)
            {
				//if( std::abs (eigvec( i2,i0 + which_band) ) > std::pow(10.,-6) )
				//if( i2 >= 10 && i2 < 31  )
				//{
				eigV.push_back( eigensolver.eigenvectors().col(i0+which_band)[i2]);
				//eigV.push_back( eigvec( i2,i0 + which_band) );
				//}
				//else eigV.push_back( std::complex<long double>(0.,0.) );

            }
            
            bloch_functions[i0].push_back( eigV );
			energies[i0].push_back( eigensolver.eigenvalues()[ i0 + which_band ] );
			//energies[i0].push_back( eigval( i0 + which_band ) );
            for(int i = 0; i < eigV.size(); i++)
            {
            //    std::cout <<i<<" "<<eigV[i]<<std::endl;
            }
			
        }
    }
}
// ========================================================================================



// ----------------------------------------------------------------------------------------
void Tbloch::print_dispersion( std::string fname ) 
{ 
	
	for(int i0 = 0; i0 < n_of_bands; ++i0)
	{
		std::fstream fs1;
		std::ostringstream ss;
		ss << i0;
		std::string si01 = ss.str();
		std::string totfn = fname +  si01 + std::string(".dat") ;
		fs1.open(totfn.c_str(), std::fstream::out); 
		fs1.precision(15);		
		
		for(int i1 = 0; i1 < energies[i0].size(); ++i1)
		{
			fs1 << i1 - ksize << " " << energies[i0][i1] << std::endl;
		}
		fs1.close();
	}

}
// ========================================================================================


    
