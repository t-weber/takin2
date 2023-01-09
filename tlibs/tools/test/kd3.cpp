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

// gcc -I. -o kd3 kd3.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include "../math/kd.h"

using namespace tl;

double fkt(double x)
{
	return std::cos(x*x);
}

int main()
{
	std::list<std::vector<double>> lst;
	for(double dx = 0.; dx<10.; dx+=0.01)
		lst.push_back(std::vector<double>({dx, fkt(dx)}));

	Kd<double> kd(lst);

	while(1)
	{
		std::string str;
		std::cout << "Point: ";
		std::getline(std::cin, str);
		std::istringstream istr(str);
		double dx, dy;
		istr >> dx >> dy;

		const std::vector<double>& vecN = kd.GetNearestNode(std::vector<double>{dx,dy});
		std::cout << "Closest Point: ";
		for(double d : vecN)
			std::cout << d << ", ";
		std::cout << std::endl;
	}

	return 0;
}
