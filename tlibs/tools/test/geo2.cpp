/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
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
