/**
 * Chrono helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date aug-2015
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

#ifndef __CHRONO_HELPERS_H__
#define __CHRONO_HELPERS_H__


#include <chrono>
#include <string>


namespace tl {

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

}
#endif
