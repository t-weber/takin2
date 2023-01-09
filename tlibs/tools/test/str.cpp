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

#include "../string/string.h"

int main()
{
	std::cout << tl::var_to_str<int>(1234567890, 10, 3) << std::endl;

	std::vector<int> vec({1,2,3,4,5,6,7,8});
	std::cout << tl::cont_to_str<decltype(vec)>(vec);
	std::cout << std::endl;

	std::cout << tl::str_is_digits(std::string("123456")) << std::endl;
	std::cout << tl::str_is_digits(std::wstring(L"123456")) << std::endl;
	std::cout << tl::str_is_digits(std::string("123a456")) << std::endl;
	std::cout << std::endl;

	std::cout << tl::ends_with<std::string>("Test123", "123") << std::endl;;
	std::cout << tl::ends_with<std::string>("Test1234", "123") << std::endl;;
	std::cout << tl::ends_with<std::string>("Test123", "21") << std::endl;;

	return 0;
}
