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

// gcc -o cofactor test/cofactor.cpp -lm -lstdc++ -std=c++11

#include <iostream>
#include "../math/linalg.h"

int main()
{
	auto m = tl::make_mat<tl::ublas::matrix<double>>({{1.,-2.,3.},{4.,5.,-6.},{-7.,8.,9.}});
	std::cout << "cof: " << tl::cofactor(m, 2,0) << std::endl;

	auto adj = tl::adjugate(m);
	std::cout << "adj: " << adj << std::endl;


	tl::ublas::matrix<double> mInv;
	tl::inverse(m, mInv);
	std::cout << "inv1: " << mInv << std::endl;

	std::cout << "inv2: " << adj/tl::determinant(m) << std::endl;

	return 0;
};
