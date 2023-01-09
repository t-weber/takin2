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

// clang -o struct test/struct.cpp -lstdc++ -lm -std=c++11

#define TLIBS_USE_GLOBAL_OPS
#include "../phys/atoms.h"
#include "../phys/mag.h"
#include <iostream>

int main()
{
	//std::cout << tl::mag_scatlen_eff<double>(4.) << std::endl;

	double a = 5.;
	double h = 0., k = 0., l = 0.;

	while(1)
	{
		std::cout << "Enter hkl: ";
		std::cin >> h >> k >> l;

		std::complex<double> F = tl::structfact(
			{
				tl::make_vec({0., 0., 0.}),
				tl::make_vec({0.5*a, 0.5*a, 0.5*a}),
			},
			tl::make_vec({h*2.*M_PI/a, k*2.*M_PI/a, l*2.*M_PI/a}),
			{1., 0.});

		std::cout << "Fn = " << F << std::endl;


		auto Fm = tl::structfact_mag(
			{
				tl::make_vec({0., 0., 0.}),
				tl::make_vec({0.5*a, 0.5*a, 0.5*a}),
			},
			{
				tl::make_vec({1., 0., 0.}),
				tl::make_vec({0., -1., 0.}),
			},
			tl::make_vec({h*2.*M_PI/a, k*2.*M_PI/a, l*2.*M_PI/a}),
			{1., 1.});

		std::cout << "Fm = " << Fm << std::endl;
	}
	return 0;
}
