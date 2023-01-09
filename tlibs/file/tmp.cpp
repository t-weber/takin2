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

#include "tmp.h"
#include "file.h"

#include <unistd.h>
#include <fcntl.h>
#include <cstdio>


namespace tl {

TmpFile::TmpFile() : m_strPrefix("tlibs_tmp"), m_iHandle(-1) {}
TmpFile::~TmpFile() { close(); }

// cygwin does not seem to have a ::mkstemp...
int TmpFile::mkstemp(std::string& strFile)
{
	std::string strRnd = tl::rand_name(s_iRndLen);
	if(!find_and_replace(strFile, std::string(s_iRndLen, 'X'), strRnd))
		return -1;

	int iFile = ::open(strFile.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
	return iFile;
}

bool TmpFile::open()
{
	namespace fs = boost::filesystem;

	fs::path pathTmpDir = fs::temp_directory_path();
	pathTmpDir /= m_strPrefix + "_" + std::string(s_iRndLen, 'X');

	m_strFile = pathTmpDir.string();
	m_iHandle = mkstemp(m_strFile);
	if(m_iHandle == -1)
		return false;

	//std::cout << "temp file: " << m_strFile << std::endl;
	return true;
}

void TmpFile::close()
{
	if(m_iHandle != -1)
		::close(m_iHandle);

	if(m_strFile != "")
		remove(m_strFile.c_str());
}

const std::string& TmpFile::GetFileName() const
{
	return m_strFile;
}

void TmpFile::SetPrefix(const char* pcStr)
{
	m_strPrefix = pcStr;
}

}
