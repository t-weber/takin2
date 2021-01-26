/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o fft_rev fft_rev.cpp -lstdc++ -std=c++11

#include "../math/fourier.h"
#include <iostream>

int main()
{
	std::cout << "1:" << std::endl;
	for(int i=0; i<1; ++i) std::cout << i << " -> " << tl::bit_reverse<int>(1,i) << std::endl;
	std::cout << std::endl;

	std::cout << "2:" << std::endl;
	for(int i=0; i<2; ++i) std::cout << i << " -> " << tl::bit_reverse<int>(2,i) << std::endl;
	std::cout << std::endl;

	std::cout << "4:" << std::endl;
	for(int i=0; i<4; ++i) std::cout << i << " -> " << tl::bit_reverse<int>(4,i) << std::endl;
	std::cout << std::endl;

	std::cout << "8:" << std::endl;
	for(int i=0; i<8; ++i) std::cout << i << " -> " << tl::bit_reverse<int>(8,i) << std::endl;
	std::cout << std::endl;

	std::cout << "16:" << std::endl;
	for(int i=0; i<16; ++i) std::cout << i << " -> " << tl::bit_reverse<int>(16,i) << std::endl;
	std::cout << std::endl;

	return 0;
}
