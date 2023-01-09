/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
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

// gcc -o tokens tokens.cpp -std=c++11 -lstdc++ -lm

#include "../string/string.h"

int main()
{
	std::string str(" \t Test<:>12<3<:>xyz    ");
	std::cout << "|" << str << "|" << std::endl;
	tl::trim(str);
	std::cout << "|" << str << "|" << std::endl;


	std::vector<std::string> vecToks;
	tl::get_tokens<std::string, std::string, std::vector<std::string>>(str, ":", vecToks);
	for(const std::string& strTok : vecToks) std::cout << strTok << std::endl;


	std::vector<std::string> vecToks2;
	tl::get_tokens_seq<std::string, std::string, std::vector>(str, "<:>", vecToks2);
	std::cout << std::endl;
	for(const std::string& strTok : vecToks2) std::cout << strTok << std::endl;



	std::string str2("TestSeparatorxyzSEPARATOR123456");

	std::vector<std::string> vecToks3;
	tl::get_tokens_seq<std::string, std::string, std::vector>(str2, "Separator", vecToks3, 0);
	std::cout << std::endl;
	for(const std::string& strTok : vecToks3) std::cout << strTok << std::endl;



	std::string str3("1234-.-5678-.-9012");

	std::vector<int> vecToks4;
	tl::get_tokens_seq<int, std::string, std::vector>(str3, "-.-", vecToks4, 0);
	std::cout << std::endl;
	for(int iTok : vecToks4) std::cout << iTok << std::endl;

	return 0;
}
