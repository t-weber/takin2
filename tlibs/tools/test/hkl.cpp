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

// gcc -o hkl hkl.cpp ../log/log.cpp -std=c++11 -lstdc++ -lm

#include "../phys/lattice.h"

using T = double;
using t_vec = tl::ublas::vector<T>;

int main()
{
	T dTh=0., dThx=0., d2Th=0., dChi=0., dPsi=0.;
	tl::Lattice<T> lat(5.,5.,5., M_PI/2., M_PI/2., M_PI/2.);

	T angle = -0.4 / 180.*M_PI;
	tl::get_euler_angles(lat,
		//tl::make_vec<t_vec>({1., 0., 0.}), tl::make_vec<t_vec>({1., 1., 0.}),
		1.4, 1., 0., 1.,
		&dTh, &dThx, &d2Th, &dChi, &dPsi,
		0., -sin(angle), cos(angle));

	tl::log_info("2theta = ", tl::r2d(d2Th));
	tl::log_info("theta = ", tl::r2d(dTh));
	tl::log_info("theta_x = ", tl::r2d(dThx));
	tl::log_info("chi = ", tl::r2d(dChi));
	tl::log_info("psi = ", tl::r2d(dPsi));

	return 0;
}
