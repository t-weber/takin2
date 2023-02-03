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
// gcc -o x3d x3d.cpp ../file/x3d.cpp -std=c++11 -lstdc++

#include "../file/x3d.h"

int main()
{
	tl::X3d x3d;

	tl::X3dTrafo *pTrafo = new tl::X3dTrafo();
	pTrafo->SetTrans(tl::make_vec({10., 0., 0.}));

	tl::X3dSphere *pSphere = new tl::X3dSphere(1.);
	pSphere->SetColor(tl::make_vec({1., 0., 0.}));
	pTrafo->AddChild(pSphere);


	tl::X3dTrafo *pTrafo2 = new tl::X3dTrafo();
	pTrafo2->SetTrans(tl::make_vec({0., 5., 0.}));

	tl::X3dSphere *pSphere2 = new tl::X3dSphere(2.);
	pSphere2->SetColor(tl::make_vec({0., 1., 0.}));
	pTrafo2->AddChild(pSphere2);

	tl::X3dPolygon *pPoly = new tl::X3dPolygon();
	pPoly->AddVertex(tl::make_vec({0., 0., 0.}));
	pPoly->AddVertex(tl::make_vec({10., 0., 0.}));
	pPoly->AddVertex(tl::make_vec({10., 10., 0.}));
	pPoly->SetColor(tl::make_vec({1., 1., 0.}));

	x3d.GetScene().AddChild(pTrafo);
	x3d.GetScene().AddChild(pTrafo2);
	x3d.GetScene().AddChild(pPoly);

	x3d.Save("sph.x3d");
	return 0;
}
