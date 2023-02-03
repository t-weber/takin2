/**
 * file helper
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 07-mar-2013
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

#ifndef __TLIB_TMPFILE__
#define __TLIB_TMPFILE__

#include "../math/rand.h"
#include "../string/string.h"

#include <boost/filesystem.hpp>


namespace tl {


// -------------------------------------------------------------------------
// simple version


template<class t_ch = char>
std::basic_string<t_ch> create_tmp_file(const char* pcPrefix)
{
	namespace fs = boost::filesystem;

	fs::path pathTmp = fs::temp_directory_path();
	pathTmp /= std::string(pcPrefix) + "_" + tl::rand_name<std::string>(8);
	std::string strTmp = pathTmp.string();

	return std::basic_string<t_ch>(strTmp.begin(), strTmp.end());
}


// -------------------------------------------------------------------------


class TmpFile
{
private:
	static const std::size_t s_iRndLen = 6;

protected:
	std::string m_strFile;
	std::string m_strPrefix;
	int m_iHandle;

public:
	TmpFile();
	~TmpFile();

	bool open();
	void close();
	const std::string& GetFileName() const;

	void SetPrefix(const char* pcStr);

	static int mkstemp(std::string& strFile);
};

}

#endif
