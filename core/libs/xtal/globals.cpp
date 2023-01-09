/**
 * globals for libcrystal
 * Warning: In Takin, the globals.* files are overriden by local ones!
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

#include "globals.h"
#include "tlibs/log/log.h"
#include "tlibs/file/file.h"
#include "tlibs/string/string.h"
#include <thread>


// -----------------------------------------------------------------------------

unsigned int g_iPrec = 6;
unsigned int g_iPrecGfx = 4;

t_real_glob g_dEps = 1e-6;
t_real_glob g_dEpsGfx = 1e-4;


// -----------------------------------------------------------------------------

unsigned int g_iMaxThreads = std::thread::hardware_concurrency();

unsigned int get_max_threads()
{
	unsigned int iMaxThreads = std::thread::hardware_concurrency();
	return std::min(iMaxThreads, g_iMaxThreads);
}

// -----------------------------------------------------------------------------


static std::vector<std::string> s_vecInstallPaths = { "." };

void add_resource_path(const std::string& strPath, bool bToBack)
{
	if(bToBack)
		s_vecInstallPaths.push_back(strPath);
	else	// insert at beginning, after "."
		s_vecInstallPaths.insert(s_vecInstallPaths.begin()+1, strPath);
}


std::string find_resource(const std::string& strFile, bool bLogErr)
{
	for(const std::string& strPrefix : s_vecInstallPaths)
	{
		std::string _strFile = strPrefix + "/" + strFile;
		//tl::log_debug("Looking for file: ", _strFile);
		if(tl::file_exists(_strFile.c_str()))
			return _strFile;
		else if(tl::file_exists((_strFile + ".gz").c_str()))
			return _strFile + ".gz";
		else if(tl::file_exists((_strFile + ".bz2").c_str()))
			return _strFile + ".bz2";
	}

	if(bLogErr)
		tl::log_err("Could not find resource file \"", strFile, "\".");
	return "";
}


std::vector<std::string> find_resource_dirs(const std::string& strDir, bool bLogErr)
{
	std::vector<std::string> vecDirs;

	for(const std::string& strPrefix : s_vecInstallPaths)
	{
		std::string _strDir = strPrefix + "/" + strDir;
		if(tl::dir_exists(_strDir.c_str()))
			vecDirs.push_back(_strDir);
	}

	if(bLogErr && vecDirs.size()==0)
		tl::log_err("Could not find resource directory \"", strDir, "\".");

	return vecDirs;
}

// -----------------------------------------------------------------------------
