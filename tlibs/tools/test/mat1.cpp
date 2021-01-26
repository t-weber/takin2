/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I. -o mat1 mat1.cpp -std=c++11 -lstdc++ -lm


#include "../math/linalg.h"
#include <iostream>

using T = double;
using t_mat = tl::ublas::matrix<T>;
using t_vec = tl::ublas::vector<T>;


int main()
{
	std::vector<t_vec> vecs =
	{
		tl::make_vec({1., 2., 3.}),
		tl::make_vec({2., 1., 4.}),
		tl::make_vec({5., 1., 3.}),
		tl::make_vec({2., 2., 3.}),
		tl::make_vec({6., 4., 2.}),
	};

	auto vecMean = tl::mean_value(vecs);
	for(t_vec& vec : vecs)
		vec -= vecMean;

	t_mat mat = tl::column_matrix(vecs);

	t_mat matCov1, matCorr1;
	std::tie(matCov1, matCorr1) = tl::covariance(vecs);
	t_mat matCov2 = tl::ublas::prod(mat, tl::ublas::trans(mat));

	std::cout << "cov(...) = " << matCov1 << std::endl;
	std::cout << "corr(...) = " << matCorr1 << std::endl;
	std::cout << "AA^t =     " <<  matCov2/matCov2(0,0)*matCov1(0,0) << std::endl;

	return 0;
}
