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

// gcc -o fs fs.cpp -lstdc++ -lboost_system -lboost_filesystem -std=c++11

#include <iostream>
#include "../file/file.h"

int main()
{
	std::cout << tl::get_file_size(std::string("test/fs.cpp")) << std::endl;

	std::vector<std::string> vec = tl::get_all_files<0>("/home/tweber/tmp");
	std::vector<std::string> vecRec = tl::get_all_files<1>("/home/tweber/tmp");

	std::cout << "\nnon-recursive: " << std::endl;
	for(const auto& file : vec)
		std::cout << file << std::endl;
	std::cout << "\n\nrecursive: " << std::endl;
	for(const auto& file : vecRec)
		std::cout << file << std::endl;
	return 0;
}
