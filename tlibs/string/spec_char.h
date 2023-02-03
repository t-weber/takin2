/**
 * Special chars
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 09-mar-14
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

#ifndef __SPEC_CHAR_H__
#define __SPEC_CHAR_H__

#include <string>
#include <unordered_map>

namespace tl {

struct SpecChar
{
	std::string strUTF8;
	std::wstring strUTF16;

	SpecChar() {}
	SpecChar(const char* pcUTF8, const wchar_t* pcUTF16)
			: strUTF8(pcUTF8), strUTF16(pcUTF16)
	{}
};

typedef std::unordered_map<std::string, SpecChar> t_mapSpecChars;


extern void init_spec_chars();
extern void deinit_spec_chars();
extern const std::string& get_spec_char_utf8(const std::string& strChar);
extern const std::wstring& get_spec_char_utf16(const std::string& strChar);

extern const t_mapSpecChars& get_spec_chars();

}

#endif
