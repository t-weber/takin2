/**
 * testing reflectivity curve
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date oct-2017
 * @license GPLv2
 *
 * gcc -std=c++14 -I../../ -o refl_tst refl_curve.cpp ../../tlibs/log/log.cpp  -lstdc++ -lboost_iostreams
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

#include "refl_curve.h"
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
	if(argc < 2)
	{
		std::cerr << "No reflectivity file given." << std::endl;
		return -1;
	}

	ReflCurve<> curv(argv[1]);
	if(!curv)
	{
		std::cerr << "Could not load file." << std::endl;
		return -1;
	}


	std::ofstream ofstr("/tmp/tst_refl.dat");
	for(double d=0; d<4.; d+=0.05)
		ofstr << d << " " << curv(d) << "\n";

	return 0;
}
