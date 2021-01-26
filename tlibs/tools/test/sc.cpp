/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
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
