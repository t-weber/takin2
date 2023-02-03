/**
 * Debug helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
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

#ifndef __TL_DEBUG_H__
#define __TL_DEBUG_H__

#include <string>
#include <boost/version.hpp>


// ----------------------------------------------------------------------------

#if BOOST_VERSION >= 105700

#include <boost/type_index.hpp>

namespace tl{

template<typename T>
std::string get_typename(bool bFull=1)
{
	boost::typeindex::type_index idx;

	if(bFull)
		idx = boost::typeindex::type_id_with_cvr<T>();
	else
		idx = boost::typeindex::type_id<T>();

	return idx.pretty_name();
}
}

#else

#include <typeinfo>

namespace tl{

template<typename T>
std::string get_typename(bool bFull=1)
{
	return std::string(typeid(T).name());
}
}

#endif

// ----------------------------------------------------------------------------


namespace tl{


#ifndef NDEBUG
	extern void log_backtrace();
#else
	static inline void log_backtrace() {}
#endif

}
#endif
