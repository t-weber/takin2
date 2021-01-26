/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o dft2 dft2.cpp ../math/fourier.cpp ../log/log.cpp -lstdc++ -lm -std=c++11

#include "../math/fourier.h"
#include <iostream>

int main()
{
	typedef std::complex<float> t_c;
	std::vector<t_c> vecIn = {t_c(1,4), t_c(2,3), t_c(3,2), t_c(4,1),
				  /*t_c(2,1), t_c(5,1), t_c(7,3), t_c(1,1)*/};

	std::cout << "in: ";
	for(const t_c& c : vecIn) std::cout << c << ", ";
	std::cout << std::endl;


	std::vector<t_c> vecOut = tl::dft_direct(vecIn, 0, 0);
	std::vector<t_c> vecOut2 = tl::fft_direct(vecIn);

	std::cout << "dft: ";
	for(const t_c& c : vecOut) std::cout << c << ", ";
	std::cout << std::endl;

	std::cout << "fft: ";
	for(const t_c& c : vecOut2) std::cout << c << ", ";
	std::cout << std::endl;

	std::cout << "ifft: ";
	std::vector<t_c> vecOut3 = tl::fft_direct(vecOut2, 1, 1);
	for(const t_c& c : vecOut3) std::cout << c << ", ";
	std::cout << std::endl;



	vecOut = tl::dft_double(vecOut);
	vecIn = tl::dft_direct(vecOut, 1, 1);

	std::cout << "idft, doubled: ";
	for(const t_c& c : vecIn) std::cout << c << ", ";
	std::cout << std::endl;

	return 0;
}
