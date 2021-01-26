/**
 * wrapper for boost r*trees
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date oct-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_RT_H__
#define __TLIBS_RT_H__

#include <vector>
#include <list>
#include <algorithm>
#include <type_traits>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>


namespace tl {

namespace geo = boost::geometry;


// ----------------------------------------------------------------------------
// static getter and setter loops for geo point
template<class t_point, class t_iter, std::size_t I, std::size_t MAX>
struct _rt_set_pt
{
	void operator()(t_point& pt, t_iter iter)
	{
		geo::set<I>(pt, *iter);

		_rt_set_pt<t_point, t_iter, I+1, MAX> s;
		s(pt, ++iter);
	}
};


template<class t_point, class t_iter, std::size_t MAX>
struct _rt_set_pt<t_point, t_iter, MAX, MAX>
{
	void operator()(t_point& pt, t_iter iter) {}
};


template<class t_point, class t_iter, std::size_t I, std::size_t MAX>
struct _rt_get_pt
{
	void operator()(const t_point& pt, t_iter iter)
	{
		*iter = geo::get<I>(pt);

		_rt_get_pt<t_point, t_iter, I+1, MAX> g;
		g(pt, ++iter);
	}
};


template<class t_point, class t_iter, std::size_t MAX>
struct _rt_get_pt<t_point, t_iter, MAX, MAX>
{
	void operator()(const t_point& pt, t_iter iter) {}
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// indexers
template<int iWhich=0, std::size_t iSize=32> struct RtIndexType {};
template<std::size_t iSize> struct RtIndexType<0, iSize>
{
	using type = geo::index::dynamic_linear;
	static type construct() { return type(iSize); }
};


template<std::size_t iSize> struct RtIndexType<1, iSize>
{
	using type = geo::index::dynamic_quadratic;
	static type construct() { return type(iSize); }
};


template<std::size_t iSize> struct RtIndexType<2, iSize>
{
	using type = geo::index::dynamic_rstar;
	static type construct() { return type(iSize); }
};


template<std::size_t iSize> struct RtIndexType<10, iSize>
{
	using type = geo::index::linear<iSize>;
	static type construct() { return type(); }
};


template<std::size_t iSize> struct RtIndexType<11, iSize>
{
	using type = geo::index::quadratic<iSize>;
	static type construct() { return type(); }
};


template<std::size_t iSize> struct RtIndexType<12, iSize>
{
	using type = geo::index::rstar<iSize>;
	static type construct() { return type(); }
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// equals that compares just the first part of a pair
template<class t_pair>
struct RtEqualsTo
{
	bool operator()(const t_pair& pair1, const t_pair& pair2) const
	{
		return geo::index::equal_to<typename t_pair::first_type>()
			(pair1.first, pair2.first);
	}
};
// ----------------------------------------------------------------------------



/**
 * R* search tree
 */
template<class T=double, std::size_t IDIM=3, std::size_t MAX_ELEMS=32, int iIdxType=2>
class Rt
{
public:
	using t_point = geo::model::point<T, IDIM, geo::cs::cartesian>;
	using t_rest = std::vector<T>;
	using t_node = std::pair<t_point, t_rest>;
	using t_idx = RtIndexType<iIdxType, MAX_ELEMS>;


protected:
	std::vector<T> m_vecMin, m_vecMax;
	using t_rt = geo::index::rtree<t_node, typename t_idx::type,
		geo::index::indexable<t_node>, RtEqualsTo<t_node>>;
	std::unique_ptr<t_rt> m_prt;


public:

	/**
	 * unloads tree
	 */
	void Unload()
	{
		m_vecMin.clear();
		m_vecMax.clear();
		m_prt.reset(nullptr);
	}

	/**
	 * creates the tree
	 */
	void Init()
	{
		Unload();
		m_prt.reset(new t_rt(t_idx::construct()));
	}


	/**
	 * insert a single point into the tree
	 */
	void InsertPoint(const std::vector<T>& vec, bool bUpdateMinMax=1)
	{
		// calc min / max
		if(bUpdateMinMax)
		{
			if(m_vecMin.size() == 0)
			{
				m_vecMin.resize(IDIM);
				m_vecMax.resize(IDIM);

				for(std::size_t i0=0; i0<IDIM; ++i0)
					m_vecMin[i0] = m_vecMax[i0] = vec[i0];
			}
			else
			{
				for(std::size_t i0=0; i0<IDIM; ++i0)
				{
					m_vecMin[i0] = std::min(m_vecMin[i0], vec[i0]);
					m_vecMax[i0] = std::max(m_vecMax[i0], vec[i0]);
				}
			}
		}


		// insert point
		t_point pt;
		_rt_set_pt<t_point, typename std::vector<T>::const_iterator, std::size_t(0), IDIM> s;
		s(pt, vec.begin());

		t_rest vecRest;
		vecRest.reserve(vec.size() - IDIM);
		std::copy(vec.begin()+IDIM, vec.end(), std::back_inserter(vecRest));

		m_prt->insert(t_node(pt, vecRest));
	}


	/**
	 * load a list of points
	 */
	void Load(const std::list<std::vector<T>>& lstPoints)
	{
		if(!m_prt) Init();

		for(const std::vector<T>& vec : lstPoints)
			InsertPoint(vec);
	}


	/**
	 * get the nodes closest to a given vector
	 */
	std::list<std::vector<T>> GetNearestNodes(const std::vector<T>& vec, std::size_t iNum=1) const
	{
		std::list<std::vector<T>> lstRet;
		if(!m_prt)
			return lstRet;

		t_point pt;
		_rt_set_pt<t_point, typename std::vector<T>::const_iterator, std::size_t(0), IDIM> s;
		s(pt, vec.begin());

		std::vector<t_node> vecRes;
		m_prt->query(geo::index::nearest(pt, iNum), std::back_inserter(vecRes));


		for(const t_node& nd : vecRes)
		{
			std::vector<T> vec;
			vec.reserve(IDIM + nd.second.size());

			_rt_get_pt<t_point, std::back_insert_iterator<std::vector<T>>, std::size_t(0), IDIM> g;
			g(nd.first, std::back_inserter(vec));

			std::copy(nd.second.begin(), nd.second.end(), std::back_inserter(vec));

			lstRet.push_back(vec);
		}

		return lstRet;
	}


	/**
	 * get the node closest to a given vector
	 */
	std::vector<T> GetNearestNode(const std::vector<T>& vec) const
	{
		std::list<std::vector<T>> lstNodes = GetNearestNodes(vec, 1);
		if(lstNodes.size() == 0)
			return std::vector<T>();
		return *lstNodes.begin();
	}


	/**
	 * tests if a vector is close to a tree node
	 */
	bool IsPointInTree(const std::vector<T>& vec,
		T eps = std::numeric_limits<T>::epsilon()) const
	{
		std::vector<T> vecNearest = GetNearestNode(vec);
		if(vecNearest.size() != IDIM)
			return false;

		for(std::size_t i=0; i<IDIM; ++i)
			if(std::abs(vec[i]-vecNearest[i]) > eps)
				return false;

		return true;
	}


	/**
	 * tests if a vector is inside the bounding box
	 */
	bool IsPointInGrid(const std::vector<T>& vec) const
	{
		for(std::size_t i=0; i<IDIM; ++i)
			if(vec[i] < m_vecMin[i] || vec[i] > m_vecMax[i])
				return false;
		return true;
	}

	Rt() { Init(); }
	Rt(std::list<std::vector<T>>& lstPoints) { Load(lstPoints); }
	~Rt() { Unload(); }
};

}
#endif
