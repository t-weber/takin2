/**
 * hashes
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date apr-2015
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

#ifndef __TLIB_HASH_H__
#define __TLIB_HASH_H__

#include <vector>
#include <algorithm>
#include <boost/functional/hash.hpp>

namespace tl {

template<typename T>
std::size_t hash(const T& t)
{
	boost::hash<T> hsh;
	return hsh(t);
}


// order of elements in container matters
template<class t_cont>
std::size_t hash_ordered(const t_cont& cont)
{
	typedef typename t_cont::const_iterator t_iter;
	typedef typename t_cont::value_type T;

	std::size_t iseed = 0;
	boost::hash<T> hsh;

	for(t_iter iter=cont.begin(); iter!=cont.end(); ++iter)
	{
		T t = *iter;
		std::size_t iHsh = hsh(t);

		boost::hash_combine(iseed, iHsh);
	}

	return iseed;
}


// order of elements in container doesn't matter
template<class t_cont>
std::size_t hash_unordered(const t_cont& cont)
{
	typedef typename t_cont::const_iterator t_iter;
	typedef typename t_cont::value_type T;

	std::vector<T> vec;

	for(t_iter iter=cont.begin(); iter!=cont.end(); ++iter)
		vec.push_back(*iter);

	std::sort(vec.begin(), vec.end());
	return hash_ordered<std::vector<T>>(vec);
}

}
#endif
