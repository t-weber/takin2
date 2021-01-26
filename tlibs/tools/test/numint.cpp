/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I. -o numint numint.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <cmath>
#include "../math/numint.h"


double fkt(double x)
{
	return x*x*x;
}


double fkt1(double x)
{
	return std::sin(x);
}

double fkt2(double x)
{
	return std::cos(x);
}

int main()
{
	std::function<double(double)> f(fkt);
	std::cout << "rect: " << tl::numint_rect(f, 3., 5., 128) << std::endl;
	std::cout << "trap: " << tl::numint_trap(f, 3., 5.) << std::endl;
	std::cout << "trapN: " << tl::numint_trapN(f, 3., 5., 128) << std::endl;
	std::cout << "simp: " << tl::numint_simp(f, 3., 5.) << std::endl;
	std::cout << "simpN: " << tl::numint_simpN(f, 3., 5., 128) << std::endl;

	std::cout << "calc: " << 0.25*5.*5.*5.*5. - 0.25*3.*3.*3.*3. << std::endl;


	std::cout << std::endl;

	for(double dX=-1.; dX<=1.; dX+=0.1)
		std::cout << "x = " << dX << ", y = "
			<< tl::convolute(std::function<double(double)>(fkt1),
				std::function<double(double)>(fkt2),
				dX, -M_PI/2., M_PI/2., 128)
			<< std::endl;

	return 0;
}
