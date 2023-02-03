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
// gcc -o geo2 geo2.cpp ../file/x3d.cpp -std=c++11 -lstdc++ -lm

#include "../file/x3d.h"
#include "../math/geo_prim.h"

int main()
{
	//tl::Dodecahedron<> prim;
	tl::TesselSphere<> prim;
	tl::X3d x3d;

	for(std::size_t iPoly=0; iPoly<prim.GetPolyCount(); ++iPoly)
	{
		tl::X3dPolygon *pPoly = new tl::X3dPolygon();
		pPoly->SetColor(tl::make_vec({1., 1., 0.}));

		for(auto vert : prim.GetPoly(iPoly))
			pPoly->AddVertex(vert);

		x3d.GetScene().AddChild(pPoly);
	}

	x3d.Save("/home/tweber/tmp/geo.x3d");
	return 0;
}
