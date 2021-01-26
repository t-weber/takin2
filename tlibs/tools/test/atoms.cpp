/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
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
