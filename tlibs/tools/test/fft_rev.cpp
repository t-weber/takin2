/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
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
