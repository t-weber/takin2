/**
 * tlibs2
 * algorithm library
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2018-2021
 * @note Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @note Forked 2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#ifndef __TLIBS2_ALGOS_H__
#define __TLIBS2_ALGOS_H__

#include <algorithm>
#include <numeric>
#include <string>
#include <chrono>
#include <vector>


namespace tl2 {

/**
 * copy algorithm with interleave
 */
template<class T1, class T2>
void copy_interleave(T1 inIter, T1 inEnd, T2 outIter, std::size_t interleave, std::size_t startskip)
{
	std::advance(inIter, startskip);

	while(std::distance(inIter, inEnd) > 0)
	{
		*outIter = *inIter;

		++outIter;
		std::advance(inIter, interleave);
	}
}


/**
 * count number of ocurrences of a sub-string in a string
 */
template<class t_str=std::string>
std::size_t count_occurrences(const t_str &str, const t_str &tok)
{
	std::size_t num = 0;
	std::size_t start = 0;
	const std::size_t len_tok = tok.length();

	while(true)
	{
		std::size_t idx = str.find(tok, start);
		if(idx == t_str::npos)
			break;

		++num;
		start += idx+len_tok;
	}

	return num;
}


/**
 * merge containers
 */
template<class t_cont>
t_cont arrayunion(const std::initializer_list<t_cont>& lst)
{
	t_cont contRet;
	for(const t_cont& cont : lst)
		contRet.insert(contRet.end(), cont.begin(), cont.end());
	return contRet;
}


/**
 * minimum of four numbers
 */
template<typename T=double>
T min4(T t1, T t2, T t3, T t4)
{
	T tmin = t1;
	tmin = std::min(tmin, t2);
	tmin = std::min(tmin, t3);
	tmin = std::min(tmin, t4);
	return tmin;
}



/**
 * like std::chrono::seconds/minutes/hours, but with variable type
 */
template<typename T = long >
using t_dur_secs = std::chrono::duration<T, std::ratio<1, 1>>;
template<typename T = long >
using t_dur_mins = std::chrono::duration<T, std::ratio<60, 1>>;
template<typename T = long >
using t_dur_hours = std::chrono::duration<T, std::ratio<60*60, 1>>;

template<typename T = long >
using t_dur_days = std::chrono::duration<T, std::ratio<60*60*24, 1>>;

template<typename T = long >
using t_dur_weeks = std::chrono::duration<T, std::ratio<60*60*24*7, 1>>;


/**
 * duration since epoch
 */
template<typename t_dur = std::chrono::seconds>
t_dur epoch_dur()
{
	namespace ch = std::chrono;
	return ch::duration_cast<t_dur>(ch::system_clock::now().time_since_epoch());
}


/**
 * seconds since epoch
 */
template<typename T=double>
T epoch()
{
	return epoch_dur<t_dur_secs<T>>().count();
}


/**
 * create a string representation of epoch
 */
template<typename T=double>
std::string epoch_to_str(T tSeconds, const char *pcFmt="%a %Y-%b-%d %H:%M:%S %Z")
{
	namespace ch = std::chrono;

	t_dur_secs<T> secs(tSeconds);
	ch::system_clock::time_point tp(ch::duration_cast<ch::seconds>(secs));

	std::time_t t = ch::system_clock::to_time_t(tp);
	std::tm tm = *std::localtime(&t);

	char cTime[256];
	std::strftime(cTime, sizeof cTime, pcFmt, &tm);
	return std::string(cTime);
}


/**
 * get the permutation of indices to sort a container
 */
template<class Comp, class t_cont = std::vector<std::size_t>>
t_cont get_perm(std::size_t num_elems, Comp comp)
{
	t_cont perm(num_elems);
	std::iota(perm.begin(), perm.end(), 0);

	std::stable_sort(perm.begin(), perm.end(), comp);
	return perm;
}


/**
 * get the permutation of indices to sort a container
 */
template<class t_cont, class t_cont_perm = std::vector<std::size_t>>
t_cont_perm get_perm(const t_cont& cont)
{
	t_cont_perm perm = get_perm(
		cont.size(),
		[&cont](std::size_t idx1, std::size_t idx2) -> bool
		{
			return cont[idx1] < cont[idx2];
		});

	return perm;
}


/**
 * reorder a vector according to a permutation
 */
template<class t_vec, class t_perm = std::vector<std::size_t>>
t_vec reorder(const t_vec& vec, const t_perm& perm)
{
	t_vec vec_new;
	vec_new.reserve(vec.size());

	for(decltype(vec.size()) i=0; i<vec.size(); ++i)
		vec_new.push_back(vec[perm[i]]);

	return vec_new;
}

}

#endif
