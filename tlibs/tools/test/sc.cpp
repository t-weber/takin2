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
// gcc -o sc sc.cpp -std=c++11 -lstdc++ -lm

#include "../math/lattice.h"
#include "../math/atoms.h"
#include "../math/nn.h"
#include "../helper/misc.h"

#include <iostream>

using namespace tl;

int main()
{
	/*std::vector<int> v{1, 2, 5, 4, 6};
	std::vector<std::size_t> vecIdx = sorted_idx(v);
	for(std::size_t idx : vecIdx)
		std::cout << idx << ", ";
	std::cout << std::endl;*/


	Lattice<double> latt(5., 5., 5., M_PI/2., M_PI/2., M_PI/2.);

	std::vector<ublas::vector<double>> vecAtoms =
	{
		make_vec({0., 0., 0.})
	};
	std::vector<std::complex<double>> vecFacts = { 1. };

	std::vector<ublas::vector<double>> vecAllAtoms;
	std::vector<std::complex<double>> vecAllFacts;
	std::tie(vecAllAtoms, vecAllFacts, std::ignore) = generate_supercell(latt, vecAtoms, vecFacts, 2);

	for(std::size_t i=0; i<vecAllAtoms.size(); ++i)
		std::cout << i << ": " << vecAllAtoms[i] << std::endl;


	std::size_t iNN=0;
	std::vector<std::vector<std::size_t>> vecIdx = get_neighbours(vecAllAtoms, vecAtoms[0]);
	for(const std::vector<std::size_t>& vecInnerIdx : vecIdx)
	{
		std::cout << "neighbour level " << iNN << ": " << std::endl;

		for(std::size_t iIdx : vecInnerIdx)
			std::cout << "\t" << vecAllAtoms[iIdx] << std::endl;

		++iNN;
	}
	return 0;
}
