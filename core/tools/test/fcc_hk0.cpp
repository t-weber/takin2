/**
 * gcc -o fcc_hk0 fcc_hk0.cpp -std=c++11 -lstdc++ -lm
 * draws 1st BZ of fcc in hk0 plane
 * @date 2015
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
	t_poly bz;
	bz.outer().push_back(t_point(-1,-1));
	bz.outer().push_back(t_point(-1,1));
	bz.outer().push_back(t_point(1,1));
	bz.outer().push_back(t_point(1,-1));
	bz.outer().push_back(t_point(-1,-1));

	t_poly bz_trans, bz1, bz2, bz3, bz4;
	geo::transform(bz, bz_trans, trafo::translate_transformer<t_real, 2, 2>(1. + 0.75*std::sqrt(2.), 0.));
	geo::transform(bz_trans, bz1, trafo::rotate_transformer<geo::degree, t_real, 2, 2>(45.));
	geo::transform(bz_trans, bz2, trafo::rotate_transformer<geo::degree, t_real, 2, 2>(-45.));
	geo::transform(bz_trans, bz3, trafo::rotate_transformer<geo::degree, t_real, 2, 2>(135.));
	geo::transform(bz_trans, bz4, trafo::rotate_transformer<geo::degree, t_real, 2, 2>(-135.));

	std::vector<t_poly> vecDiff1, vecDiff2, vecDiff3, vecDiff4;
	geo::difference(bz,bz1, vecDiff1);
	geo::difference(vecDiff1[0],bz2, vecDiff2);
	geo::difference(vecDiff2[0],bz3, vecDiff3);
	geo::difference(vecDiff3[0],bz4, vecDiff4);

	geo::svg_mapper<t_point> svg(std::cout, 1000, 1000);

/*	svg.add(bz); svg.map(bz, "stroke:rgb(0,0,0); fill:rgb(0,0,255)");
	svg.add(bz1); svg.map(bz1, "stroke:rgb(0,0,0); fill:rgb(0,255,0)");
	svg.add(bz2); svg.map(bz2, "stroke:rgb(0,0,0); fill:rgb(0,255,0)");
	svg.add(bz3); svg.map(bz3, "stroke:rgb(0,0,0); fill:rgb(0,255,0)");
	svg.add(bz4); svg.map(bz4, "stroke:rgb(0,0,0); fill:rgb(0,255,0)");*/

	for(const t_poly& poly : vecDiff4)
	{
		//geo::transform(poly, poly, trafo::scale_transformer<t_real, 2,2>(1e-5, 1e-5));
		svg.add(poly);
		svg.map(poly, "stroke:rgb(0,0,0); fill:rgb(0,0,255)");
	}



	std::vector<t_point> vecPts =
	{
		t_point(0, 0),		// Gamma
		t_point(0.75, 0.75),	// K
		t_point(0.75, -0.75),	// K
		t_point(-0.75, 0.75),	// K
		t_point(-0.75, -0.75),	// K

		t_point(0., 1.),	// X
		t_point(0., -1.),	// X
		t_point(1., 0.),	// X
		t_point(-1., 0.),	// X
	};

	for(const t_point& pt : vecPts)
	{
		svg.add(pt);
		svg.map(pt, "stroke:rgb(0,0,0); fill:rgb(0,0,0)", 10);
	}

	return 0;
}
