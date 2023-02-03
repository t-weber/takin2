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

// gcc -I. -o numint numint.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <cmath>
#include "../math/numint.h"


double fkt(double x)
{
	return x*x*x;
}


double fkt1(double x)
{
	return std::sin(x);
}

double fkt2(double x)
{
	return std::cos(x);
}

int main()
{
	std::function<double(double)> f(fkt);
	std::cout << "rect: " << tl::numint_rect(f, 3., 5., 128) << std::endl;
	std::cout << "trap: " << tl::numint_trap(f, 3., 5.) << std::endl;
	std::cout << "trapN: " << tl::numint_trapN(f, 3., 5., 128) << std::endl;
	std::cout << "simp: " << tl::numint_simp(f, 3., 5.) << std::endl;
	std::cout << "simpN: " << tl::numint_simpN(f, 3., 5., 128) << std::endl;

	std::cout << "calc: " << 0.25*5.*5.*5.*5. - 0.25*3.*3.*3.*3. << std::endl;


	std::cout << std::endl;

	for(double dX=-1.; dX<=1.; dX+=0.1)
		std::cout << "x = " << dX << ", y = "
			<< tl::convolute(std::function<double(double)>(fkt1),
				std::function<double(double)>(fkt2),
				dX, -M_PI/2., M_PI/2., 128)
			<< std::endl;

	return 0;
}
