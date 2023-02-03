/**
 * gcc -o fcc_hhl fcc_hhl.cpp -std=c++11 -lstdc++ -lm
 * draws 1st BZ of fcc in hhl plane
 * @date 2016
 * @license GPLv2
 * @author Tobias Weber <tobias.weber@tum.de>
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

#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>
namespace geo = boost::geometry;
namespace trafo = geo::strategy::transform;

using t_real = double;
using t_point = geo::model::point<t_real, 2, geo::cs::cartesian>;
using t_poly = geo::model::polygon<t_point>;

int main()
{
	t_real s2 = std::sqrt(2.);

	t_poly bz;
	bz.outer().push_back(t_point(-0.25*s2,	-1.));
	bz.outer().push_back(t_point(0.25*s2,	-1.));
	bz.outer().push_back(t_point(0.75*s2,	0.));
	bz.outer().push_back(t_point(0.25*s2,	1.));
	bz.outer().push_back(t_point(-0.25*s2,	1.));
	bz.outer().push_back(t_point(-0.75*s2,	0.));
	bz.outer().push_back(t_point(-0.25*s2,	-1.));


	geo::svg_mapper<t_point> svg(std::cout, 1000, 1000);

	svg.add(bz);
	svg.map(bz, "stroke:rgb(0,0,0); fill:rgb(0,0,255)");


	std::vector<t_point> vecPts =
	{
		t_point(0, 0),		// Gamma

		t_point(0.75*s2, 0.),	// K
		t_point(-0.75*s2, 0.),	// K

		t_point(0., 1.),		// X
		t_point(0., -1.),		// X

		t_point(0.5*s2, 0.5),		// L
		t_point(-0.5*s2, 0.5),		// L
		t_point(0.5*s2, -0.5),		// L
		t_point(-0.5*s2, -0.5),		// L
	};

	for(const t_point& pt : vecPts)
	{
		svg.add(pt);
		svg.map(pt, "stroke:rgb(0,0,0); fill:rgb(0,0,0)", 10);
	}

	return 0;
}
