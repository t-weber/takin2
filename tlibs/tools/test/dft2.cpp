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
