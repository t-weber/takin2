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

// gcc -o rt rt.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <vector>
#include <list>
#include "../math/rt.h"

int main()
{
	std::list<std::vector<int>> lst =
	{
		{1,9},
		{2,8},
		{3,7},
		{4,6},
		{5,5},
		{6,4},
		{7,3},
		{8,2},
		{9,1},
		{10,0},
		{20,15}
	};

	tl::Rt<int, 2> rt(lst);

	const std::vector<int>& vecN = rt.GetNearestNode(std::vector<int>{20,20});
	std::cout << "nearest: ";
	for(int i : vecN)
		std::cout << i << ", ";
	std::cout << std::endl;
}
