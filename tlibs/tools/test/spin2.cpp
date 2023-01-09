/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 * clang -I/usr/include/lapacke -o spin2 ../test/spin2.cpp ../log/log.cpp -lstdc++ -std=c++11 -lm
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

#include "../phys/spin.h"
#include <iostream>

using t_real = double;
using namespace tl;

int main()
{
	t_real J, j1, j2, m1, m2;
	std::cout << "J = "; std::cin >> J;
	std::cout << "j1 = "; std::cin >> j1;
	std::cout << "j2 = "; std::cin >> j2;
	std::cout << "m1 = "; std::cin >> m1;
	std::cout << "m2 = "; std::cin >> m2;

	std::cout << CG_coeff(J, j1, j2, m1, m2) << std::endl;
	return 0;
}
