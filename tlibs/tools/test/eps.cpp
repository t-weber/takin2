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

// gcc -o eps eps.cpp -lstdc++ -std=c++11

#include "../math/math.h"
#include "../math/linalg.h"
#include <complex>
#include <iostream>

int main()
{
	std::cout << tl::get_epsilon<double>() << std::endl;
	std::cout << tl::get_epsilon<std::complex<double>>() << std::endl;
	std::cout << tl::get_epsilon<tl::ublas::matrix<std::complex<double>>>() << std::endl;

	std::cout << tl::float_equal<double>(0., 1.) << std::endl;
	std::cout << tl::float_equal<std::complex<double>>(
		std::complex<double>(0.,1.), std::complex<double>(1.,1.))
		<< std::endl;
	return 0;
}
