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

// gcc -o atoms atoms.cpp -std=c++11 -lstdc++ -lm

#include "../math/linalg.h"
#include "../math/atoms.h"
#include <vector>
#include <iostream>

using namespace tl;

int main()
{
	typedef ublas::matrix<double> t_mat;
	typedef ublas::vector<double> t_vec;

	std::vector<t_mat> vecMat =
	{
		make_mat({{1.,0.},{0.,1.}}),
		make_mat({{0.,1.},{1.,0.}}),
		make_mat({{0.,1.},{1.,0.}}),
	};

	t_vec vecAtom = make_vec({1., 2.});

	std::vector<t_vec> vecResult =
		generate_atoms<t_mat, t_vec, std::vector>(vecMat, vecAtom);


	for(const t_vec& vecRes : vecResult)
		std::cout << vecRes << std::endl;
	return 0;
}
