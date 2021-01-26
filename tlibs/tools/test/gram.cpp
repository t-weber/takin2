/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// clang -o gram gram.cpp -std=c++11 -lstdc++ -lm

#include "../math/linalg.h"

namespace ublas = boost::numeric::ublas;

int main()
{
	ublas::vector<double> vec1 = tl::make_vec({1.,2.,3.});
	ublas::vector<double> vec2 = tl::make_vec({2.,-3.,4.});
	ublas::vector<double> vec3 = tl::make_vec({3.,4.,5.});

	std::vector<ublas::vector<double>> vecsIn = {vec1, vec2, vec3};
	std::vector<ublas::vector<double>> vecsOut = tl::gram_schmidt(vecsIn);

	for(const ublas::vector<double>& vec : vecsOut)
		std::cout << vec << std::endl;

	std::cout << "v0 * v1 = " << ublas::inner_prod(vecsOut[0], vecsOut[1]) << std::endl;
	std::cout << "v0 * v2 = " << ublas::inner_prod(vecsOut[0], vecsOut[2]) << std::endl;
	std::cout << "v1 * v2 = " << ublas::inner_prod(vecsOut[1], vecsOut[2]) << std::endl;



	// --------------------------------------------------------------------
	std::cout << "\n";

	ublas::vector<double> vecA = tl::make_vec({-2.,0.5,2.});
	ublas::vector<double> vecB = tl::make_vec({1.,-2.,1.});

	std::vector<ublas::vector<double>> vecsOutC = tl::gram_schmidt({vecA, vecB});

	for(const ublas::vector<double>& vec : vecsOutC)
		std::cout << vec << ", ";
	std::cout << std::endl;


        ublas::vector<double> vecUp = tl::cross_3(vecA, vecB);
        vecB = tl::cross_3(vecUp, vecA);

        vecA /= ublas::norm_2(vecA);
        vecB /= ublas::norm_2(vecB);
        vecUp /= ublas::norm_2(vecUp);

	for(const ublas::vector<double>& vec : {vecA, vecB, vecUp})
		std::cout << vec << ", ";
	std::cout << std::endl;

	return 0;
}
