/**
 * array helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
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

#ifndef __TLIBS_ARRAYS_H__
#define __TLIBS_ARRAYS_H__

#include <cstddef>
#include <vector>
#include <list>
#include <array>
#include <tuple>
#include <algorithm>


namespace tl {


/**
 * Minimalistic wrapper for plain old arrays
 */
template<class T>
class wrapper_array
{
	public:
		using value_type = T;
		using size_type = std::size_t;

		using iterator = T*;
		using const_iterator = const T*;

		using reference = T&;
		using const_reference = const T&;


	protected:
		T* m_pt = nullptr;
		size_type m_len = 0;

	public:
		wrapper_array(T* p, size_type len)
			: m_pt(p), m_len(len)
		{}

		wrapper_array() = default;
		~wrapper_array() = default;

		size_type size() const { return m_len; }

		iterator begin() { return m_pt; }
		iterator end() { return m_pt+m_len; }
		const_iterator begin() const { return m_pt; }
		const_iterator end() const { return m_pt+m_len; }

		reference operator[](size_type i) { return m_pt[i]; }
		const_reference operator[](size_type i) const { return m_pt[i]; }
};



// ----------------------------------------------------------------------------
// conversions

/**
 * Converts a t_cont_in to a t_cont_out
 */
template<template<class...> class t_cont_out, class t_cont_in>
t_cont_out<typename t_cont_in::value_type> convert_containers(const t_cont_in& cont1)
{
	using T = typename t_cont_in::value_type;

	t_cont_out<T> cont2;

	for(const T& val : cont1)
		cont2.push_back(val);

	return cont2;
}


template<typename T>
std::list<T> vector_to_list(const std::vector<T>& vec)
{
	return convert_containers<std::list, std::vector<T>>(vec);
}


template<typename T>
T* vec_to_array(const std::vector<T>& vec)
{
	T* t_arr = new T[vec.size()];

	std::size_t i=0;
	for(const T& t : vec)
		t_arr[i++] = t;

	return t_arr;
}

// ----------------------------------------------------------------------------



/**
 * sort tuple-vector
 */
template<const std::size_t isortidx,
	class... Ts,
	template<class...> class t_cont = std::vector>
void sorttuples(t_cont<std::tuple<Ts...> >& vec)
{
	std::sort(vec.begin(), vec.end(),
		[](const std::tuple<Ts...>& tup1, const std::tuple<Ts...>& tup2) -> bool
		{ return std::get<isortidx>(tup1) < std::get<isortidx>(tup2);});
}


// ----------------------------------------------------------------------------

template<class t_cont>
t_cont arrayunion(const std::initializer_list<t_cont>& lst)
{
	t_cont contRet;
	for(const t_cont& cont : lst)
		contRet.insert(contRet.end(), cont.begin(), cont.end());
	return contRet;
}


// ----------------------------------------------------------------------------


/**
 * converts container t_from to container t_to
 */
template<class t_to, class t_from, template<class...> class t_cont=std::vector,
	bool bIsEqu = std::is_same<t_from, t_to>::value>
struct container_cast
{
	using t_from_al = std::allocator<t_from>;
	using t_to_al = std::allocator<t_to>;

	t_cont<t_to, t_to_al> operator()(const t_cont<t_from, t_from_al>& vec) const
	{
		t_cont<t_to, t_to_al> vecTo;
		for(const t_from& t : vec)
			vecTo.push_back(t_to(t));
		return vecTo;
	}
};

/**
 * if t_from == t_to  ->  return reference
 */
template<class t_to, class t_from, template<class...> class t_cont>
struct container_cast<t_to, t_from, t_cont, 1>
{
	using t_from_al = std::allocator<t_from>;
	using t_to_al = std::allocator<t_to>;

	const t_cont<t_to, t_to_al>& operator()(const t_cont<t_from, t_from_al>& vec) const
	{ return vec; }

	t_cont<t_to, t_to_al>& operator()(t_cont<t_from, t_from_al>& vec) const
	{ return vec; }
};


// ----------------------------------------------------------------------------


/**
 * converts container of container with t_from to container of container with t_to
 */
template<class t_to, class t_from, template<class...> class t_cont=std::vector,
	bool bIsEqu = std::is_same<t_from, t_to>::value>
struct container2_cast
{
	using t_from_al = std::allocator<t_from>;
	using t_to_al = std::allocator<t_to>;
	using t_from_al2 = std::allocator<t_cont<t_from, t_from_al>>;
	using t_to_al2 = std::allocator<t_cont<t_to, t_to_al>>;

	t_cont<t_cont<t_to, t_to_al>, t_to_al2>
	operator()(const t_cont<t_cont<t_from, t_from_al>, t_from_al2>& vec) const
	{
		t_cont<t_cont<t_to, t_to_al>, t_to_al2> vecvecTo;
		for(const auto& vecInner : vec)
		{
			t_cont<t_to, t_to_al> vecTo;
			for(const t_from& t : vecInner)
				vecTo.push_back(t_to(t));

			vecvecTo.emplace_back(std::move(vecTo));
		}

		return vecvecTo;
	}
};


/**
 * if t_from == t_to  ->  return reference
 */
template<class t_to, class t_from, template<class...> class t_cont>
struct container2_cast<t_to, t_from, t_cont, 1>
{
	using t_from_al = std::allocator<t_from>;
	using t_to_al = std::allocator<t_to>;
	using t_from_al2 = std::allocator<t_cont<t_from, t_from_al>>;
	using t_to_al2 = std::allocator<t_cont<t_to, t_to_al>>;

	const t_cont<t_cont<t_to, t_to_al>, t_to_al2>&
	operator()(const t_cont<t_cont<t_from, t_from_al>, t_from_al2>& vec) const
	{ return vec; }

	t_cont<t_cont<t_to, t_to_al>, t_to_al2>&
	operator()(t_cont<t_cont<t_from, t_from_al>, t_from_al2>& vec) const
	{ return vec; }
};


// ----------------------------------------------------------------------------


template<std::size_t iWhichElem = 0,
	template<class...> class t_cont = std::vector,
	class t_pairvec = std::vector<std::pair<double,double>>>
struct vec_from_pairvec
{
	using t_pair = typename t_pairvec::value_type;
	using t_val = typename std::tuple_element<iWhichElem, t_pair>::type;

	t_cont<t_val> operator()(const t_pairvec& vec) const
	{
		t_cont<t_val> vecVal;
		vecVal.reserve(vec.size());

		for(const t_pair& elem : vec)
			vecVal.push_back(std::get<iWhichElem>(elem));

		return vecVal;
	}
};

// ----------------------------------------------------------------------------


}
#endif
