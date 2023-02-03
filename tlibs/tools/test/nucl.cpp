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

// gcc -I. -o nucl nucl.cpp -std=c++11 -lstdc++ -lm

#include "../math/nuclear.h"
#include <fstream>

int main()
{
	std::ofstream ofstr("collisions.dat");

	for(double A=1.001; A<6.; A+=0.001)
	{
		tl::t_energy_si<double> E_from = 2. * tl::get_one_MeV<double>();
		tl::t_energy_si<double> E_to = 25. * tl::get_one_meV<double>();

		ofstr << std::setw(25) << std::left << A << " "
			<< tl::mean_collisions(A, E_from, E_to) << std::endl;
	}

	return 0;
}
