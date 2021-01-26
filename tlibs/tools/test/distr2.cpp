/**
 * tlibs test file
 * Test distributions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date apr-17
 * @license GPLv2 or GPLv3
 */

// gcc -std=c++11 -o distr2 distr2.cpp -lstdc++ -lm


#include "../math/stat.h"
#include <iostream>
#include <vector>

using t_real = double;


int main()
{
	std::vector<t_real> v({1, 5, 3, 2, 7, 9, 15, 50});

	t_real dDof = v.size()-1;
	tl::t_student_dist<t_real> t(dDof);
	std::cout << "t(p=0.9) = " << t.cdf_inv(0.9) << std::endl;

	t_real dMean, dStd, dConf;
	std::tie(dMean, dStd, dConf) = tl::confidence<t_real, decltype(v)>(v, 0.8);
	std::cout << "80% confidence: <v> = " << dMean << " +- " << dConf << std::endl;

	std::tie(dMean, dStd, dConf) = tl::confidence<t_real, decltype(v)>(v, 0.9);
	std::cout << "90% confidence: <v> = " << dMean << " +- " << dConf << std::endl;

	std::tie(dMean, dStd, dConf) = tl::confidence<t_real, decltype(v)>(v, 0.99);
	std::cout << "99% confidence: <v> = " << dMean << " +- " << dConf << std::endl;

	std::cout << "S = " << tl::entropy(t, -10., 10., 100) << std::endl;


	std::cout << "\n";

	tl::t_binomial_dist<t_real> b(50, 0.5);
	std::cout << "S = " << tl::entropy(b, 50) << std::endl;;

	tl::t_normal_dist<t_real> n(0., 1.);
	std::cout << "S = " << tl::entropy(n, -10., 10., 100) << std::endl;;


	std::cout << "\n";
	std::cout << "ordered, repeating particles: " << tl::combinatorics<t_real>(10, 5, 1, 1) << std::endl;
	std::cout << "ordered, non-repeating particles: " << tl::combinatorics<t_real>(10, 5, 1, 0) << std::endl;
	std::cout << "non-ordered, repeating particles: " << tl::combinatorics<t_real>(10, 5, 0, 1) << std::endl;
	std::cout << "non-ordered, non-repeating particles: " << tl::combinatorics<t_real>(10, 5, 0, 0) << std::endl;

	std::cout << "\n";
	std::cout << "Boltzons: " << tl::boltzons_in_niveaus<t_real>(2,2) << std::endl;
	std::cout << "Bosons: " << tl::bosons_in_niveaus<t_real>(2,2) << std::endl;
	std::cout << "Fermions: " << tl::fermions_in_niveaus<t_real>(2,2) << std::endl;

	return 0;
}
