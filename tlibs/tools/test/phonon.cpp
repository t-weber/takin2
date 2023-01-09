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

// gcc -I. -o phonon test/phonon.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <fstream>
#include "../math/neutrons.h"


int main()
{
	std::ofstream ofstr("phonontst.gpl");

	double dT = 120.;
	double dE0 = 2.;
	double dHWHM = 1.;
	double dAmp = 1.;
	double dOffs = 0.;

	ofstr << "plot \"-\" pointtype 7\n";
	for(double dE=-5.; dE<5.; dE+=0.1)
	{
		double dS = tl::DHO_model<double>(dE, dT, dE0, dHWHM, dAmp, dOffs);
		ofstr << dE << "\t" << dS << "\n";
	}
	ofstr << "end\n";

	return 0;
}
