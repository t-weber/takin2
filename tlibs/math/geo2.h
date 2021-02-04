/**
 * advanced geometry helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 25-apr-2020
 * @license GPLv2 or GPLv3
 *
 * @see https://github.com/t-weber/misc/blob/master/geo/qhulltst.cpp
 *
 * References:
 *   - http://www.qhull.org/html/qh-code.htm#cpp
 *   - https://github.com/qhull/qhull/tree/master/src/libqhullcpp
 *   - https://github.com/qhull/qhull/blob/master/src/qhulltest/Qhull_test.cpp
 */

#ifndef __TLIBS_GEO2_H__
#define __TLIBS_GEO2_H__

#include "geo.h"

#ifndef NO_QHULL
	#include <memory>

	#include <Qhull.h>
	#include <QhullFacetList.h>
	#include <QhullVertexSet.h>
#endif


namespace tl {

namespace ublas = boost::numeric::ublas;
namespace math = boost::math;



#ifndef NO_QHULL

using t_real_qh = double;

template<class t_vec = ublas::vector<double>,
template<class...> class t_cont = std::vector,
class T = typename t_vec::value_type>
t_cont<t_cont<t_vec>> get_convexhull_qh(const t_cont<t_vec>& vecVerts)
{
	t_cont<t_cont<t_vec>> vecPolys;
	const t_vec vecCentre = mean_value(vecVerts);

	// copy vertices
	int dim = vecVerts[0].size();
	std::size_t len = vecVerts.size()*dim;
	std::unique_ptr<t_real_qh[]> mem{new t_real_qh[len]};

	std::size_t i=0;
	for(const t_vec& vert : vecVerts)
	{
		for(int d=0; d<dim; ++d)
		{
			mem[i] = t_real_qh(vert[d]);
			++i;
		}
	}


	orgQhull::Qhull qhull{"tlibs_geo2", dim, int(vecVerts.size()), mem.get(), ""};

	orgQhull::QhullFacetList facets = qhull.facetList();
	//std::cout << "Convex hull has " << facets.size() << " facets.";

	for(orgQhull::QhullLinkedList<orgQhull::QhullFacet>::iterator iter=facets.begin();
		iter!=facets.end(); ++iter)
	{
		// triangulable?
		if(iter->isUpperDelaunay())
			continue;

		t_cont<t_vec> vecPoly;
		orgQhull::QhullVertexSet vertices = iter->vertices();
		for(orgQhull::QhullSet<orgQhull::QhullVertex>::iterator iterVertex=vertices.begin();
			iterVertex!=vertices.end(); ++iterVertex)
		{
			orgQhull::QhullPoint point = (*iterVertex).point();
			t_vec vecPoint(dim);
			for(int i=0; i<dim; ++i)
				vecPoint[i] = T(point[i]);

			vecPoly.emplace_back(std::move(vecPoint));
		}

		sort_poly_verts<t_vec, t_cont, T>(vecPoly, vecCentre);
		vecPolys.emplace_back(std::move(vecPoly));
	}

	// too few polygons => remove polyhedron
	if(vecPolys.size() < 3)
		vecPolys = decltype(vecPolys){};

	return vecPolys;
}

#endif



/**
 * creates polyhedron faces out of vertices
 */
template<class t_vec = ublas::vector<double>,
template<class...> class t_cont = std::vector,
class T = typename t_vec::value_type>
t_cont<t_cont<t_vec>> get_convexhull(
	const t_cont<t_vec>& vecVerts, T eps = tl::get_epsilon<T>())
{
#ifdef NO_QHULL
	// note that his function is more restrictive than a convex hull:
	// it searches for polyhedra and does not permit additional vertices
	// inside the polyhedron
	return verts_to_polyhedron<t_vec, t_cont, T>(vecVerts, eps);
#else
	return get_convexhull_qh<t_vec, t_cont, T>(vecVerts);
#endif
}


}

#endif
