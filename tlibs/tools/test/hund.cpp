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

// clang -o hund hund.cpp -std=c++11 -lstdc++

#include <iostream>
#include "../math/term.h"

int main()
{
	//auto tup = tl::get_orbital<std::string>("2p3");
	//std::cout << std::get<0>(tup) << ", " << std::get<1>(tup) << ", " << std::get<2>(tup) << std::endl;


	double S,L,J;
	std::uint16_t l, es;

	std::cout << "l = "; std::cin >> l;
	std::cout << "Number of electrons: "; std::cin >> es;

	//std::tie(S,L,J) = tl::hund(std::string("1s2 2p3"));
	std::tie(S,L,J) = tl::hund(l, es);
	std::string strTerm = tl::get_termsymbol(S,L,J);

	std::cout << strTerm << "\n";
	std::cout << "S = " << S << "\n";
	std::cout << "L = " << L << "\n";
	std::cout << "J = " << J << "\n";
	std::cout.flush();

	return 0;
}
