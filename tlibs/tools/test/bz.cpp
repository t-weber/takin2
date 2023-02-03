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

// gcc -DNO_LAPACK -o bz bz.cpp ../log/log.cpp ../file/x3d.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <fstream>
#include "../phys/bz.h"
#include "../file/x3d.h"

using T = double;
using t_vec = tl::Brillouin3D<T>::t_vec<T>;

int main()
{
	tl::Brillouin3D<T> bz;

	std::ifstream ifstr("/home/tweber/tmp/bz.dat");
	t_vec vecCentr(3);
	ifstr >> vecCentr[0] >> vecCentr[1] >> vecCentr[2];

	bz.SetCentralReflex(vecCentr);

	while(1)
	{
		t_vec vec(3);
		ifstr >> vec[0] >> vec[1] >> vec[2];
		if(!ifstr)
			break;

		bz.AddReflex(vec);
	}

	bz.CalcBZ();
	std::cout << bz.GetVertices().size() << " vertices." << std::endl;
	std::cout << bz.GetPolys().size() << " polygons." << std::endl;


	tl::X3d x3d;
	for(const t_vec& vec : bz.GetVertices())
	{
		tl::X3dTrafo *pTrafo = new tl::X3dTrafo();
		pTrafo->SetTrans(vec - bz.GetCentralReflex());

		tl::X3dSphere *pSphere = new tl::X3dSphere(0.025);
		pSphere->SetColor(tl::make_vec({1., 0., 0.}));
		pTrafo->AddChild(pSphere);

		x3d.GetScene().AddChild(pTrafo);
	}

	for(const std::vector<t_vec>& vecPoly : bz.GetPolys())
	{
		//tl::X3dPolygon *pPoly = new tl::X3dPolygon();
		tl::X3dLines *pPoly = new tl::X3dLines();
		pPoly->SetColor(tl::make_vec({1., 1., 0.}));

		for(const t_vec& vec : vecPoly)
			pPoly->AddVertex(vec - bz.GetCentralReflex());

		x3d.GetScene().AddChild(pPoly);
	}



	tl::Plane<T> plane(
		bz.GetCentralReflex(),  	    // vec0
		tl::make_vec({1., 0., 0.}));	// norm
	auto tupLinesandVerts = bz.GetIntersection(plane);
	std::cout << "slice: " << std::get<0>(tupLinesandVerts).size() << " lines." << std::endl;
	std::cout << "slice: " << std::get<1>(tupLinesandVerts).size() << " vertices." << std::endl;

	for(const t_vec& vec : std::get<1>(tupLinesandVerts))
	{
		tl::X3dTrafo *pTrafo = new tl::X3dTrafo();
		pTrafo->SetTrans(vec - bz.GetCentralReflex());

		tl::X3dSphere *pSphere = new tl::X3dSphere(0.05);
		pSphere->SetColor(tl::make_vec({0., 0., 1.}));
		pTrafo->AddChild(pSphere);

		x3d.GetScene().AddChild(pTrafo);
	}


	x3d.Save("/home/tweber/tmp/bz.x3d");	

	return 0;
}
