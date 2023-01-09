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

// gcc -I../.. -o tst_hkl tst_hkl.cpp ../../tlibs/log/log.cpp -lstdc++ -std=c++11 -lm 

#include <iostream>
#include "tlibs/math/lattice.h"

using namespace tl;

typedef ublas::vector<double> t_vec;


int main()
{
	double a = 5.;
	double b = 4.;
	double c = 3.;
	double alpha = 90.;
	double beta = 90.;
	double gamma = 60.;

	double dKi = 2.66;
	double dKf = 2.66;

	double dh = 2.;
	double dk = 0.;
	double dl = 1.;

	t_vec vec1 = make_vec({1., -1., 0.});
	t_vec vec2 = make_vec({1.,  1., 1.});

	Lattice<double> lat(a,b,c, alpha/180.*M_PI,beta/180.*M_PI,gamma/180.*M_PI);
	//Lattice<double> rec = lat.GetRecip();
	//Lattice<double> rec_aligned = rec.GetAligned();

	std::cout << "lattice: " << lat << std::endl;
	//std::cout << "recip: " << rec << std::endl;
	//std::cout << "recip (aligned): " << rec_aligned << std::endl;

	std::cout << "vec1 = " << vec1 << std::endl;
	std::cout << "vec2 = " << vec2 << std::endl;


	double dTheta, d2Theta;

	try
	{
		get_tas_angles<double>(lat, vec1,vec2, dKi,dKf, dh,dk,dl, 1, &dTheta, &d2Theta);
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
		return -1;
	}

	std::cout << "theta = " << dTheta/M_PI*180. << " deg" << std::endl;
	std::cout << "2theta = " << d2Theta/M_PI*180. << " deg" << std::endl;
	return 0;
}
