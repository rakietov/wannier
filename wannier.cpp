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
Twannier::Twannier( Tbloch bl1  ) 
{
    ksize = bl1.return_ksize();
    lsize = bl1.return_lsize();
    xx_step = bl1.return_xx_step();
    x_max = bl1.return_x_max();
    n_of_bands = bl1.return_n_of_bands();  
    lattice = bl1.return_lattice();    
    bloch_functions = bl1.return_bloch_functions();

// setup XX matrix:

    //std::cout << "setup XX matrix" << std::endl;
    Eigen::Matrix<  std::complex<long double>, Eigen::Dynamic, Eigen::Dynamic> XXm;
    XXm.resize( 2*ksize + 1, 2*ksize + 1);



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
                        std::complex<long double> tmp =  std::complex<long double>(double(std::pow((-1.), i1-i2+n1-n2))/(( 
                       double(i1-i2)/(2.*ksize+1.)+ n1-n2 )));
                        tmp *= std::complex<long double>(0.,0.5)  * std::conj( bloch_functions[i2+ksize][n2] ) *  bloch_functions[i1+ksize][n1];
                        mat_ele += tmp;
                    }
                }
            }
            XXm(i1+ksize,i2+ksize) = mat_ele;
            //std::cout << i1<<" "<<i2<<" "<<mat_ele << std::endl;
        }
    }
    //std::cout << "setup XX matrix DONE" << std::endl;

	Eigen::SelfAdjointEigenSolver< Eigen::Matrix<std::complex<long double>, Eigen::Dynamic, Eigen::Dynamic> > eigensolver(XXm);


    std::fstream fs1, fs2;
    fs1.open("tmp.dat", std::fstream::out);	
    std::cout <<eigensolver.eigenvalues();
    fs1 << eigensolver.eigenvectors().col( ksize );
    fs1.close();
            
    fs2.open("tmp.dat", std::fstream::in);	

    std::vector< std::complex<long double> > eigV;        
    for(int i2 = 0; i2 < 2*ksize+5; i2++)
    {
        std::string str1;
        fs2 >> str1;
        std::istringstream isstmp(str1); 
        std::complex<long double> tmp(0.0,0.0);
        bool flag = (isstmp >> tmp);
        if(flag) eigV.push_back(tmp);
    }
    fs2.close();
    //std:: cout <<std::endl<< eigV.size()<<std:: endl;

    for(int ix = 0; ix < 10000; ix++)
    {
        double x = ix/100. - 50.;
        std::complex<long double> wannier_x (0., 0.);
        for(int i1= -ksize; i1 < ksize +1; ++i1)
        {
            std::complex<long double> tmp = 0.;
            for(int i2 = -lsize; i2 < lsize+1; ++i2)
            {
                tmp += bloch_functions[i1+ksize][i2+lsize] * std::exp( std::complex<long double>(0., (i2 + i1 *1./(2.*ksize+1.))* x)   );
                //std::cout << tmp << std::endl;
            }
            wannier_x += eigV[i1+ksize] * tmp;
        }
        std::cout << x << " " << wannier_x.real()<<" "<<wannier_x.imag()<< std::endl;
    }


}
// ========================================================================================

