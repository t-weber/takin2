/**
 * tlibs
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2012-2020
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

#include "version.h"
#include <string>

namespace tl {


const char* get_tlibs_version()
{
	return TLIBS_VERSION;
}

const char* get_tlibs_infos()
{
	return "This is the tlibs, a physical-mathematical C++ template library.\n"
		"Written by Tobias Weber <tobias.weber@tum.de>, 2012 - 2021.\n"
		"License: GPLv2 or GPLv3.";
}

/**
 * Check if supplied string matches with compiled one
 * (to check if .so library and header files match)
 */
bool check_tlibs_version(const char* pcHdrVer)
{
	return std::string(TLIBS_VERSION) == std::string(pcHdrVer);
}


}
