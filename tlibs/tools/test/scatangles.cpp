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

// gcc -o scatangles test/scatangles.cpp -lm -lstdc++ -std=c++11

#include <iostream>
#include "../math/neutrons.h"

int main()
{
	using T = float;

	T dKi = 1.4;
	T dKf = 1.4;
	T dQ = 2.5;

	while(1)
	{
		std::cout << "ki = "; std::cin >> dKi;
		std::cout << "kf = "; std::cin >> dKf;
		std::cout << "Q = "; std::cin >> dQ;

		tl::t_wavenumber_si<T> ki = dKi / tl::get_one_angstrom<T>();
		tl::t_wavenumber_si<T> kf = dKf / tl::get_one_angstrom<T>();
		tl::t_wavenumber_si<T> Q = dQ / tl::get_one_angstrom<T>();

		T kiQ = tl::get_angle_ki_Q(ki, kf, Q, 1, 0) / tl::get_one_radian<T>();
		T kfQ = tl::get_angle_kf_Q(ki, kf, Q, 1, 0) / tl::get_one_radian<T>();

		std::cout << "kiQ = " << tl::r2d(kiQ) << std::endl;
		std::cout << "kfQ = " << tl::r2d(kfQ) << std::endl;

		std::cout << std::endl;
	}

	return 0;
}
