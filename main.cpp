#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <complex>

#include "lattice.h"
#include "bloch.h"
#include "wannier.h"

#include <unistd.h>

using namespace std;

int main(){


		double a = 0.0 ;    
		double V0 = 40.;
		auto phase = std::exp(std::complex<long double>(0.,0.5*3.14159265));
		std::vector< std::complex<long double> > coeffs = { 
        std::complex<long double>(V0*((1.+a)/2.),0.), 
        std::complex<long double>((-a/4.)*V0 ,0.)*phase  , 
        std::complex<long double>( -(1./4.)*V0,0.)*phase*phase  };
    
	    Tlattice lat1( coeffs );
		Tbloch blo1(lat1, 40, 20, 100., 2, 1);


		Twannier wan1( blo1,0);


/*
	fstream fs1, fs2;
	fs1.open("s_energies.dat",std::fstream::out);

	fs2.open("p_energies.dat",std::fstream::out);

	for(int i = 0; i < 100; ++i)
	{
		double a = 0.01 * i;    
		double V0 = 40.;
		auto phase = std::exp(std::complex<long double>(0.,0.5*3.14159265));
		std::vector< std::complex<long double> > coeffs = { 
        std::complex<long double>(V0*((1.+a)/2.),0.), 
        std::complex<long double>((-a/4.)*V0 ,0.)*phase  , 
        std::complex<long double>( -(1./4.)*V0,0.)*phase*phase  };
    
	    Tlattice lat1( coeffs );
		Tbloch blo1(lat1, 80, 40, 100., 2, 0);
		//blo1.print_dispersion("dispersion");

		std::ostringstream ss1;
		ss1 << a;
		string s_a = ss1.str();

		fstream fs1, fs2;


		Twannier wan1( blo1,0);
		wan1.print_wannier( (string("wan_s0_a_") + s_a + string(".dat")).c_str() );
	
		Twannier wan2( blo1,1);
		wan2.print_wannier( (string("wan_s1_a_") + s_a + string(".dat")).c_str() );

		fs1 << a << " " << wan1.calc_energy() << " " << wan2.calc_energy() <<std:: endl;

		Tbloch blo2(lat1, 80, 40, 100., 2, 1);
		//blo1.print_dispersion("dispersion");

		Twannier wan3( blo2,0);
		wan3.print_wannier( (string("wan_p0_a_") + s_a + string(".dat")).c_str() );
	
		Twannier wan4( blo2,1);
		wan4.print_wannier( (string("wan_p1_a_") + s_a + string(".dat")).c_str() );

		fs2 << a << " " << wan3.calc_energy() << " " << wan4.calc_energy() <<std:: endl;


	}
*/
	return 0;
}
















