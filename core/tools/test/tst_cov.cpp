/**
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

// clang -I../.. -I. -o tst_cov tst_cov.cpp -lstdc++ -lm -std=c++11

#include <iostream>
#include <vector>
#include <initializer_list>
#include <type_traits>
#include "tlibs/math/linalg.h"

using namespace tl;

int main()
{
	std::vector<ublas::vector<double>> vals
	{
		make_vec({1.0, 2.0, 3.0}),
		make_vec({1.1, 2.2, 3.4}),
		make_vec({0.9, 1.8, 2.9}),
		make_vec({1.0, 2.1, 3.1}),
		make_vec({1.4, 2.1, 3.6}),
	};

	std::vector<double> vecProb = { 1./5., 1./5., 1./5., 1./5., 1./5.};

	std::cout << "Measured values:\n";
	std::copy(vals.begin(), vals.end(), std::ostream_iterator<decltype(vals[0])>(std::cout,"\n"));
	std::cout << std::endl;

	ublas::matrix<double> mat = covariance(vals, &vecProb);
	std::cout << "Covariance matrix:\n" << mat << std::endl;

	return 0;
}
