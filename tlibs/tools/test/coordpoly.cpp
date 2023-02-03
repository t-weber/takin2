/**
 * tests calculation of coordination polyhedra
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

// g++ -std=c++11 -I. -I/usr/include/libqhullcpp -DNO_LAPACK -o coordpoly tools/test/coordpoly.cpp log/log.cpp file/x3d.cpp -lqhull_r -lqhullcpp

#include <iostream>
#include <fstream>
#include "math/geo2.h"
#include "file/x3d.h"

using T = double;
using t_vec = tl::ublas::vector<T>;


void SaveX3d(const std::vector<std::vector<t_vec>>& vecPolys,
	const std::vector<t_vec>& vecVertices, const char* pcFile)
{
	t_vec vecCentral = tl::zero_v<t_vec>(3);

	tl::X3d<T> x3d;
	for(const t_vec& vec : vecVertices)
	{
		tl::X3dTrafo<T> *pTrafo = new tl::X3dTrafo<T>();
		pTrafo->SetTrans(vec - vecCentral);

		tl::X3dSphere<T> *pSphere = new tl::X3dSphere<T>(0.025);
		pSphere->SetColor(tl::make_vec({1., 0., 0.}));
		pTrafo->AddChild(pSphere);

		x3d.GetScene().AddChild(pTrafo);
	}

	for(const auto& vecPoly : vecPolys)
	{
		tl::X3dPolygon<T> *pPoly = new tl::X3dPolygon<T>();
		//tl::X3dLines *pPoly = new tl::X3dLines();
		pPoly->SetColor(tl::make_vec({1., 1., 0.}));

		for(const t_vec& vec : vecPoly)
			pPoly->AddVertex(vec - vecCentral);

		x3d.GetScene().AddChild(pPoly);
	}

	x3d.Save(pcFile);
}


int main()
{
	std::ifstream ifstr("/tmp/coord.dat");
	std::vector<t_vec> vecVertices;

	while(1)
	{
		t_vec vec(3);
		ifstr >> vec[0] >> vec[1] >> vec[2];
		if(!ifstr)
			break;

		vecVertices.push_back(vec);
	}

	std::vector<std::vector<t_vec>> vecPolys = tl::verts_to_polyhedron<t_vec, std::vector, T>(vecVertices, 1e-4);
	std::vector<std::vector<t_vec>> vecPolys_qh = tl::get_convexhull_qh<t_vec, std::vector, T>(vecVertices);

	SaveX3d(vecPolys, vecVertices, "/tmp/coord.x3d");
	SaveX3d(vecPolys_qh, vecVertices, "/tmp/coord_qh.x3d");

	return 0;
}
