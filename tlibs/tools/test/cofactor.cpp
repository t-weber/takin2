/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o cofactor test/cofactor.cpp -lm -lstdc++ -std=c++11

#include <iostream>
#include "../math/linalg.h"

int main()
{
	auto m = tl::make_mat<tl::ublas::matrix<double>>({{1.,-2.,3.},{4.,5.,-6.},{-7.,8.,9.}});
	std::cout << "cof: " << tl::cofactor(m, 2,0) << std::endl;

	auto adj = tl::adjugate(m);
	std::cout << "adj: " << adj << std::endl;


	tl::ublas::matrix<double> mInv;
	tl::inverse(m, mInv);
	std::cout << "inv1: " << mInv << std::endl;

	std::cout << "inv2: " << adj/tl::determinant(m) << std::endl;

	return 0;
};
