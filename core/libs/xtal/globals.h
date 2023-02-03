/**
 * globals
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 20-mar-2015
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#ifndef __XTL_GLOBALS_H__
#define __XTL_GLOBALS_H__

#include <vector>
#include <string>


//using t_real_glob = float;
using t_real_glob = double;


extern unsigned int g_iPrec;
extern unsigned int g_iPrecGfx;

extern t_real_glob g_dEps;
extern t_real_glob g_dEpsGfx;

extern unsigned int g_iMaxThreads;

extern void add_resource_path(const std::string& strPath, bool bToBack=1);
extern std::string find_resource(const std::string& strFile, bool bLogErr=1);
extern std::vector<std::string> find_resource_dirs(const std::string& strDir, bool bLogErr=1);

extern unsigned int get_max_threads();


#endif
