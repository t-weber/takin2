/**
 * tlibs2
 * bit manipulations
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2012-2021
 * @license GPLv3, see 'LICENSE' file
 * @note Forked on 4-Sep-2021 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
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

#ifndef __TLIBS2_BITS_H__
#define __TLIBS2_BITS_H__


namespace tl2 {

/*
 * count how many bits are needed for the given number
 */
template<typename T = std::size_t>
T count_needed_bits(T imax)
{
	T inum = 0;

	for(; imax!=0; imax>>=1)
		++inum;

	return inum;
}


/**
 * reverse the bits of the given number
 */
template<typename T = std::size_t>
T bit_reverse(T imax, T inum)
{
	if(imax<2) return inum;

	T irev = 0;
	T ibitcnt = count_needed_bits<T>(imax)-2;

	for(T i=1; i<imax; i<<=1)
	{
		if(inum & i)
			irev |= (1 << ibitcnt);
		--ibitcnt;
	}

	return irev;
}


template<typename T = std::size_t,
	template<class...> class t_cont = std::vector>
t_cont<T> bit_reverse_indices(T imax)
{
	t_cont<T> vec;
	vec.reserve(imax);

	for(T i=0; i<imax; ++i)
		vec.push_back(bit_reverse<T>(imax, i));

	return vec;
}

}

#endif
