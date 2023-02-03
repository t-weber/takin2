/**
 * Exceptions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 04-mar-2014
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

#ifndef __MY_EXCEPT_H__
#define __MY_EXCEPT_H__

#include <exception>
#include <string>

namespace tl {

class Err : public std::exception
{
	protected:
		std::string m_strErr;

	public:
		Err(const std::string& strErr, bool bErrStr=0) noexcept
			: m_strErr((bErrStr? "Exception: " : "") + strErr)
		{}

		Err(const char* pcErr) noexcept : Err(std::string(pcErr))
		{}

		virtual ~Err() noexcept
		{}

		virtual const char* what() const noexcept override
		{
			return m_strErr.c_str();
		}
};

}
#endif
