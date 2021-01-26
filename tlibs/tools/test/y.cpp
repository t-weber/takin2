/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */
// clang -o y y.cpp -std=c++11 -lstdc++ -lm

#include "../math/math.h"
#include <fstream>
#include <iostream>
#include <iomanip>

int main()
{
	const int N = 64;
	int l=0, m=0;

	std::cout << "l = "; std::cin >> l;
	std::cout << "m = "; std::cin >> m;

	std::ofstream ofstr("y.dat");
	ofstr << "# spherical harmonic with l=" << l << " and m=" << m << "\n";
	ofstr << "#" << std::setw(15) << "theta" << std::setw(16) << "phi";
	ofstr << std::setw(16) << "real" << std::setw(16) << "imag" << "\n";

	for(double th=0; th<M_PI; th+=M_PI/double(N))
	for(double ph=0; ph<2.*M_PI; ph+=2.*M_PI/double(N))
	{
		std::complex<double> c = tl::Ylm<double>(l,m,th,ph);

		ofstr << std::setw(16) << th << std::setw(16) << ph;
		ofstr << std::setw(16) << c.real() << std::setw(16) << c.imag();
		ofstr << "\n";
	}

	return 0;
}
