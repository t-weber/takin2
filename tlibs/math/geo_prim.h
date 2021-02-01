/**
 * geometric primitives
 * @author: Tobias Weber <tobias.weber@tum.de>
 * @date: 14-may-2017
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_GEO_PRIM_H__
#define __TLIBS_GEO_PRIM_H__

#include "geo.h"
#include "stat.h"


namespace tl {

namespace ublas = boost::numeric::ublas;
namespace math = boost::math;


template<class t_vec = ublas::vector<double>>
class GeometricPrimitive
{
public:
	using t_vertices = std::vector<t_vec>;
	using t_polyindices = std::vector<std::vector<std::size_t>>;
	using t_polyindex = typename t_polyindices::value_type;

protected:
	virtual t_vertices& GetVertices() = 0;
	virtual t_polyindices& GetPolyIndices() = 0;

public:
	virtual std::size_t GetVertexCount() const = 0;
	virtual std::size_t GetPolyCount() const = 0;

	virtual const t_vec& GetVertex(std::size_t iVert) const = 0;
	virtual const t_polyindex& GetPolyIndex(std::size_t iPoly) const = 0;


	/**
	 * polygons comprising the solid
	 */
	virtual t_vertices GetPoly(std::size_t iPoly) const
	{
		std::vector<t_vec> vecPoly;

		for(std::size_t iIdx : GetPolyIndex(iPoly))
			vecPoly.push_back(GetVertex(iIdx));

		return vecPoly;
	}


	/**
	 * normal vectors to polygon faces
	 */
	virtual t_vec GetPolyNormal(std::size_t iPoly) const
	{
		t_vec vecNorm;

		t_polyindex vecPolyIdx = GetPolyIndex(iPoly);
		if(vecPolyIdx.size() < 3)	// too few vertices
			return vecNorm;

		t_vec vec1 = GetVertex(vecPolyIdx[1]) - GetVertex(vecPolyIdx[0]);
		t_vec vec2 = GetVertex(vecPolyIdx[2]) - GetVertex(vecPolyIdx[1]);
		vecNorm = cross_3(vec1, vec2);
		vecNorm /= veclen(vecNorm);

		return vecNorm;
	}


	/**
	 * tesselate using a central position
	 */
	virtual void SubdividePolysInMiddle()
	{
		t_vertices& vecVertices = GetVertices();
		t_polyindices& vecvecPolyIndices = GetPolyIndices();
		t_polyindices vecvecNewPolyIndices;

		// iterate over all polys
		for(const auto& vecPolyIndices : vecvecPolyIndices)
		{
			std::size_t iNumVerts = vecPolyIndices.size();
			std::vector<t_vec> vecVerts;

			for(std::size_t iVert=0; iVert<iNumVerts; ++iVert)
				vecVerts.push_back(GetVertex(vecPolyIndices[iVert]));

			// add a new vertex in the middle of the polygon
			t_vec vecVertMid = mean_value(vecVerts);
			vecVertices.push_back(vecVertMid);
			std::size_t iIdxVertMid = vecVertices.size()-1;

			// new vertex indices
			for(std::size_t iVert=0; iVert<iNumVerts; ++iVert)
			{
				std::size_t iVertIdx = vecPolyIndices[iVert];
				std::size_t iNextVertIdx = vecPolyIndices[(iVert+1) % iNumVerts];

				std::vector<std::size_t> vecNewVerts = { iVertIdx, iNextVertIdx, iIdxVertMid };
				vecvecNewPolyIndices.emplace_back(std::move(vecNewVerts));
			}
		}

		// replace poly indices
		vecvecPolyIndices = std::move(vecvecNewPolyIndices);
	}


	virtual void SubdividePolysInMiddle(std::size_t iIters)
	{
		for(std::size_t iIter=0; iIter<iIters; ++iIter)
			SubdividePolysInMiddle();
	}


	/**
	 * tesselate along polygon edges
	 */
	virtual void SubdividePolysAlongEdges()
	{
		t_vertices& vecVertices = GetVertices();
		t_polyindices& vecvecPolyIndices = GetPolyIndices();
		t_polyindices vecvecNewPolyIndices;

		// iterate over all polys
		for(const auto& vecPolyIndices : vecvecPolyIndices)
		{
			std::size_t iNumVerts = vecPolyIndices.size();
			if(iNumVerts != 3)	// only for triangles!
				break;

			t_vec vecMidEdge1 =
				(GetVertex(vecPolyIndices[0]) + GetVertex(vecPolyIndices[1])) / 2.;
			t_vec vecMidEdge2 =
				(GetVertex(vecPolyIndices[1]) + GetVertex(vecPolyIndices[2])) / 2.;
			t_vec vecMidEdge3 =
				(GetVertex(vecPolyIndices[2]) + GetVertex(vecPolyIndices[0])) / 2.;

			// add a new vertices along the edges
			vecVertices.push_back(vecMidEdge1);
			std::size_t iIdxMidEdge1 = vecVertices.size()-1;
			vecVertices.push_back(vecMidEdge2);
			std::size_t iIdxMidEdge2 = vecVertices.size()-1;
			vecVertices.push_back(vecMidEdge3);
			std::size_t iIdxMidEdge3 = vecVertices.size()-1;

			// new polygon indices
			vecvecNewPolyIndices.emplace_back(std::vector<std::size_t>
				({ vecPolyIndices[0], iIdxMidEdge1, iIdxMidEdge3 }));
			vecvecNewPolyIndices.emplace_back(std::vector<std::size_t>
				({ iIdxMidEdge1, vecPolyIndices[1], iIdxMidEdge2 }));
			vecvecNewPolyIndices.emplace_back(std::vector<std::size_t>
				({ iIdxMidEdge3, iIdxMidEdge2, vecPolyIndices[2] }));
			vecvecNewPolyIndices.emplace_back(std::vector<std::size_t>
				({ iIdxMidEdge1, iIdxMidEdge2, iIdxMidEdge3 }));
		}

		// replace poly indices
		vecvecPolyIndices = std::move(vecvecNewPolyIndices);
	}


	virtual void SubdividePolysAlongEdges(std::size_t iIters)
	{
		for(std::size_t iIter=0; iIter<iIters; ++iIter)
			SubdividePolysAlongEdges();
	}

};


// ----------------------------------------------------------------------------


/**
 * Tetrahedron
 * @see e.g.: https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec = ublas::vector<double>>
class Tetrahedron : public GeometricPrimitive<t_vec>
{
public:
	using t_vertices = typename GeometricPrimitive<t_vec>::t_vertices;
	using t_polyindices = typename GeometricPrimitive<t_vec>::t_polyindices;
	using t_polyindex = typename GeometricPrimitive<t_vec>::t_polyindex;

protected:
	// vertices
	t_vertices m_vecVertices =
	{
		make_vec<t_vec>({  1,  1,  1 }),	// 0
		make_vec<t_vec>({  1, -1, -1 }),	// 1
		make_vec<t_vec>({ -1, -1,  1 }),	// 2
		make_vec<t_vec>({ -1,  1, -1 }),	// 3
	};

	// polygons
	t_polyindices m_vecPolyIndices =
	{
		{ 0, 2, 1 }, { 0, 1, 3 },
		{ 0, 3, 2 }, { 1, 2, 3 },
	};

protected:
	virtual t_vertices& GetVertices() override { return m_vecVertices; }
	virtual t_polyindices& GetPolyIndices() override { return m_vecPolyIndices; }

public:
	Tetrahedron() = default;
	~Tetrahedron() = default;

	virtual std::size_t GetVertexCount() const override
	{ return m_vecVertices.size(); }

	virtual std::size_t GetPolyCount() const override
	{ return m_vecPolyIndices.size(); }

	virtual const t_vec& GetVertex(std::size_t iVert) const override
	{ return m_vecVertices[iVert]; }

	virtual const t_polyindex& GetPolyIndex(std::size_t iPoly) const override
	{ return m_vecPolyIndices[iPoly]; }
};


// ----------------------------------------------------------------------------


/**
 * Cube
 * @see e.g.: https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec = ublas::vector<double>>
class Cube : public GeometricPrimitive<t_vec>
{
public:
	using t_vertices = typename GeometricPrimitive<t_vec>::t_vertices;
	using t_polyindices = typename GeometricPrimitive<t_vec>::t_polyindices;
	using t_polyindex = typename GeometricPrimitive<t_vec>::t_polyindex;

protected:
	// vertices
	t_vertices m_vecVertices =
	{
		make_vec<t_vec>({  1,  1,  1 }),	// 0
		make_vec<t_vec>({  1,  1, -1 }),	// 1
		make_vec<t_vec>({  1, -1, -1 }),	// 2
		make_vec<t_vec>({  1, -1,  1 }),	// 3
		make_vec<t_vec>({ -1,  1,  1 }),	// 4
		make_vec<t_vec>({ -1,  1, -1 }),	// 5
		make_vec<t_vec>({ -1, -1, -1 }),	// 6
		make_vec<t_vec>({ -1, -1,  1 }),	// 7
	};

	// polygons
	t_polyindices m_vecPolyIndices =
	{
		{ 4, 5, 6, 7 } /*-x*/, { 7, 6, 2, 3 } /*-y*/, { 1, 2, 6, 5 } /*-z*/,
		{ 3, 2, 1, 0 } /*+x*/, { 5, 4, 0, 1 } /*+y*/, { 4, 7, 3, 0 } /*+z*/,
	};

protected:
	virtual t_vertices& GetVertices() override { return m_vecVertices; }
	virtual t_polyindices& GetPolyIndices() override { return m_vecPolyIndices; }

public:
	Cube() = default;
	~Cube() = default;

	virtual std::size_t GetVertexCount() const override
	{ return m_vecVertices.size(); }

	virtual std::size_t GetPolyCount() const override
	{ return m_vecPolyIndices.size(); }

	virtual const t_vec& GetVertex(std::size_t iVert) const override
	{ return m_vecVertices[iVert]; }

	virtual const t_polyindex& GetPolyIndex(std::size_t iPoly) const override
	{ return m_vecPolyIndices[iPoly]; }
};


// ----------------------------------------------------------------------------


/**
 * Octahedron
 * @see e.g.: https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec = ublas::vector<double>>
class Octahedron : public GeometricPrimitive<t_vec>
{
public:
	using t_vertices = typename GeometricPrimitive<t_vec>::t_vertices;
	using t_polyindices = typename GeometricPrimitive<t_vec>::t_polyindices;
	using t_polyindex = typename GeometricPrimitive<t_vec>::t_polyindex;

protected:
	// vertices
	t_vertices m_vecVertices =
	{
		make_vec<t_vec>({  1,  0,  0 }),	// 0
		make_vec<t_vec>({ -1,  0,  0 }),	// 1
		make_vec<t_vec>({  0,  1,  0 }),	// 2
		make_vec<t_vec>({  0, -1,  0 }),	// 3
		make_vec<t_vec>({  0,  0,  1 }),	// 4
		make_vec<t_vec>({  0,  0, -1 }),	// 5
	};

	// polygons
	t_polyindices m_vecPolyIndices =
	{
		{ 0, 2, 4 }, { 0, 5, 2 }, { 0, 4, 3 }, { 0, 3, 5 },
		{ 1, 4, 2 }, { 1, 2, 5 }, { 1, 3, 4 }, { 1, 5, 3 },
	};

protected:
	virtual t_vertices& GetVertices() override { return m_vecVertices; }
	virtual t_polyindices& GetPolyIndices() override { return m_vecPolyIndices; }

public:
	Octahedron() = default;
	~Octahedron() = default;

	virtual std::size_t GetVertexCount() const override
	{ return m_vecVertices.size(); }

	virtual std::size_t GetPolyCount() const override
	{ return m_vecPolyIndices.size(); }

	virtual const t_vec& GetVertex(std::size_t iVert) const override
	{ return m_vecVertices[iVert]; }

	virtual const t_polyindex& GetPolyIndex(std::size_t iPoly) const override
	{ return m_vecPolyIndices[iPoly]; }
};


// ----------------------------------------------------------------------------


/**
 * Icosahedron
 * @see e.g.: https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec = ublas::vector<double>>
class Icosahedron : public GeometricPrimitive<t_vec>
{
public:
	using t_vertices = typename GeometricPrimitive<t_vec>::t_vertices;
	using t_polyindices = typename GeometricPrimitive<t_vec>::t_polyindices;
	using t_polyindex = typename GeometricPrimitive<t_vec>::t_polyindex;

protected:
	// Golden Ratio
	/*static constexpr*/ const typename t_vec::value_type s_g = 0.5 + 0.5*std::sqrt(5.);

	// vertices
	t_vertices m_vecVertices =
	{
		make_vec<t_vec>({  0,  1,  s_g }),	// 0
		make_vec<t_vec>({  0,  1, -s_g }),	// 1
		make_vec<t_vec>({  0, -1,  s_g }),	// 2
		make_vec<t_vec>({  0, -1, -s_g }),	// 3

		make_vec<t_vec>({  1,  s_g,  0 }),	// 4
		make_vec<t_vec>({  1, -s_g,  0 }),	// 5
		make_vec<t_vec>({ -1,  s_g,  0 }),	// 6
		make_vec<t_vec>({ -1, -s_g,  0 }),	// 7

		make_vec<t_vec>({  s_g,  0,  1 }),	// 8
		make_vec<t_vec>({  s_g,  0, -1 }),	// 9
		make_vec<t_vec>({ -s_g,  0,  1 }),	// 10
		make_vec<t_vec>({ -s_g,  0, -1 }),	// 11
	};

	// polygons
	t_polyindices m_vecPolyIndices =
	{
		{ 0, 10, 2 }, { 0, 2, 8 }, { 0, 8, 4 }, { 0,  4, 6 }, { 0,  6, 10 },	// upper cap
		{ 3,  5, 7 }, { 3, 9, 5 }, { 3, 1, 9 }, { 3, 11, 1 }, { 3,  7, 11 },	// lower cap
		{ 10, 7, 2 }, { 2, 5, 8 }, { 8, 9, 4 }, { 4,  1, 6 }, { 6, 11, 10 },	// sides
		{  7, 5, 2 }, { 5, 9, 8 }, { 9, 1, 4 }, { 1, 11, 6 }, { 11, 7, 10 },	// sides
	};

protected:
	virtual t_vertices& GetVertices() override { return m_vecVertices; }
	virtual t_polyindices& GetPolyIndices() override { return m_vecPolyIndices; }

public:
	Icosahedron() = default;
	~Icosahedron() = default;

	virtual std::size_t GetVertexCount() const override
	{ return m_vecVertices.size(); }

	virtual std::size_t GetPolyCount() const override
	{ return m_vecPolyIndices.size(); }

	virtual const t_vec& GetVertex(std::size_t iVert) const override
	{ return m_vecVertices[iVert]; }

	virtual const t_polyindex& GetPolyIndex(std::size_t iPoly) const override
	{ return m_vecPolyIndices[iPoly]; }
};


// ----------------------------------------------------------------------------


/**
 * Dodecahedron
 * @see e.g.: https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec = ublas::vector<double>>
class Dodecahedron : public GeometricPrimitive<t_vec>
{
public:
	using t_vertices = typename GeometricPrimitive<t_vec>::t_vertices;
	using t_polyindices = typename GeometricPrimitive<t_vec>::t_polyindices;
	using t_polyindex = typename GeometricPrimitive<t_vec>::t_polyindex;

protected:
	// Golden Ratio
	/*static constexpr*/ const typename t_vec::value_type s_g = 0.5 + 0.5*std::sqrt(5.);

	// vertices
	t_vertices m_vecVertices =
	{
		make_vec<t_vec>({  0,  1./s_g,  s_g }),	// 0
		make_vec<t_vec>({  0,  1./s_g, -s_g }),	// 1
		make_vec<t_vec>({  0, -1./s_g,  s_g }),	// 2
		make_vec<t_vec>({  0, -1./s_g, -s_g }),	// 3

		make_vec<t_vec>({  1./s_g,  s_g,  0 }),	// 4
		make_vec<t_vec>({  1./s_g, -s_g,  0 }),	// 5
		make_vec<t_vec>({ -1./s_g,  s_g,  0 }),	// 6
		make_vec<t_vec>({ -1./s_g, -s_g,  0 }),	// 7

		make_vec<t_vec>({  s_g,  0,  1./s_g }),	// 8
		make_vec<t_vec>({  s_g,  0, -1./s_g }),	// 9
		make_vec<t_vec>({ -s_g,  0,  1./s_g }),	// 10
		make_vec<t_vec>({ -s_g,  0, -1./s_g }),	// 11

		make_vec<t_vec>({  1,  1,  1 }),	// 12
		make_vec<t_vec>({  1,  1, -1 }),	// 13
		make_vec<t_vec>({  1, -1, -1 }),	// 14
		make_vec<t_vec>({  1, -1,  1 }),	// 15

		make_vec<t_vec>({ -1,  1,  1 }),	// 16
		make_vec<t_vec>({ -1,  1, -1 }),	// 17
		make_vec<t_vec>({ -1, -1, -1 }),	// 18
		make_vec<t_vec>({ -1, -1,  1 }),	// 19
	};

	// polygons
	t_polyindices m_vecPolyIndices =
	{
		{ 16, 10, 19,  2, 0 }, { 19, 7,  5, 15,  2 }, { 15, 8, 12, 0, 2 }, { 12, 4, 6, 16, 0 },	// top cap
		{ 18, 11, 17,  1, 3 }, {  6, 4, 13,  1, 17 }, { 13, 9, 14, 3, 1 }, { 14, 5, 7, 18, 3 },	// bottom cap
		{ 19, 10, 11, 18, 7 }, { 16, 6, 17, 11, 10 }, { 15, 5, 14, 9, 8 }, { 12, 8, 9, 13, 4 },	// sides
	};

protected:
	virtual t_vertices& GetVertices() override { return m_vecVertices; }
	virtual t_polyindices& GetPolyIndices() override { return m_vecPolyIndices; }

public:
	Dodecahedron() = default;
	~Dodecahedron() = default;

	virtual std::size_t GetVertexCount() const override
	{ return m_vecVertices.size(); }

	virtual std::size_t GetPolyCount() const override
	{ return m_vecPolyIndices.size(); }

	virtual const t_vec& GetVertex(std::size_t iVert) const override
	{ return m_vecVertices[iVert]; }

	virtual const t_polyindex& GetPolyIndex(std::size_t iPoly) const override
	{ return m_vecPolyIndices[iPoly]; }
};


// ----------------------------------------------------------------------------


/**
 * tessellated sphere
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_underlying_solid = Icosahedron>
class TesselSphere : public t_underlying_solid<t_vec>
{
public:
	TesselSphere(typename t_vec::value_type dRad = 1., std::size_t iSubdivisions=2)
	{
		t_underlying_solid<t_vec>::SubdividePolysAlongEdges(iSubdivisions);

		for(t_vec& vecVertex : t_underlying_solid<t_vec>::m_vecVertices)
		{
			vecVertex /= veclen(vecVertex);
			vecVertex *= dRad;
		}
	}

	~TesselSphere() = default;
};

}
#endif
