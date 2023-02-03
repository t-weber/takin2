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

// gcc -I. -o fermi fermi.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <fstream>
#include "../math/neutrons.h"

using T = double;


int main()
{
	std::ofstream ofstr("fermitst.dat");
	T dMu = 5.;

	for(T dT=0.; dT<20.; dT+=5.)
	{
		for(T dE=0.; dE<10.; dE+=0.1)
		{
			T n = tl::fermi<T>(dE, dMu, dT);
			ofstr << dE << "\t" << n << "\n";
		}
		ofstr << "e\n";
	}

	return 0;
}
