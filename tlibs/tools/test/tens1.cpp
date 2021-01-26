 /**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I. -o tens1 tens1.cpp ../log/log.cpp -std=c++11 -lstdc++ -lm


#include "../math/tensor.h"
#include <iostream>

using T = double;
using t_mat = tl::ublas::matrix<T>;
using t_vec = tl::ublas::vector<T>;


int main()
{
	T dAngle = 35.;
	T s = std::sin(dAngle/180.*M_PI);
	T c = std::cos(dAngle/180.*M_PI);

	std::vector<t_vec> vecsBaseCov =
	{
		tl::make_vec({1., 0., 1.}) / std::sqrt(2.),
		tl::make_vec({ c,  s, 0.}),
		tl::make_vec({1., 1., -2.}) / std::sqrt(6.),
	};

	t_mat matGContra;
	t_mat matGCov = tl::make_metric_cov(vecsBaseCov);
	if(!tl::inverse(matGCov, matGContra))
		return -1;

	std::cout << "g_cov = " << matGCov << std::endl;
	std::cout << "g_contra = " << matGContra << std::endl;

	std::cout << "det(g_cov) = " << tl::determinant(matGCov) << std::endl;
	std::cout << "det(g_contra) = " << tl::determinant(matGContra) << std::endl;
	std::cout << std::endl;

	
	{
		t_mat matBaseCov = tl::column_matrix(vecsBaseCov);
		//t_mat matBaseContra = tl::column_matrix(vecsBaseContra);
		t_mat matBaseContra = tl::ublas::prod(matBaseCov, matGContra);
		std::cout << "base_cov = " << matBaseCov << std::endl;
		std::cout << "base_contra = " << matBaseContra << std::endl;
		std::cout << "det(base_cov) = " << tl::determinant(matBaseCov) << std::endl;
		std::cout << "det(base_contra) = " << tl::determinant(matBaseContra) << std::endl;

		// check cov[i]*contra[j] = delta[i,j]
		std::cout << "base cov * base contra = " << tl::ublas::prod(tl::ublas::trans(matBaseContra), matBaseCov) << std::endl;
		std::cout << std::endl;
	}


	{
		t_vec vec0 = tl::make_vec({1.,0.,0.});
		t_vec vec1 = tl::make_vec({0.,1.,0.});

		t_vec vec0cov = tl::ublas::prod(matGCov, vec0);
		t_vec vec1cov = tl::ublas::prod(matGCov, vec1);

		std::cout << "v1 * v2 = " << tl::inner_prod(matGCov, vec0, vec1) << std::endl;
		std::cout << "v1 * v2 = " << tl::inner_prod(matGContra, vec0cov, vec1cov) << std::endl;

		std::cout << "angle(v1,v2) = " << tl::vec_angle(matGCov, vec0, vec1)/M_PI*180. << std::endl;

		t_vec vecCross = tl::cross_prod_contra(matGCov, vec0, vec1);
		std::cout << "v1 x v2 (contra) = " << vecCross << std::endl;
		std::cout << "v1 x v2 (cov) = " << tl::ublas::prod(matGCov, vecCross) << std::endl;
		std::cout << "v1 x v2 (coord) = " << tl::cross_3(vecsBaseCov[0], vecsBaseCov[1]) << std::endl;
		// sign depends on det(base)
		std::cout << "coord test = " << vecCross[0]*vecsBaseCov[0] + vecCross[1]*vecsBaseCov[1] + vecCross[2]*vecsBaseCov[2] << std::endl;
	}

	{
		std::cout << std::endl;
		auto pairPerm = tl::count_permutations({3,1,2}, {1,2,3});
		std::cout << "permutation = " << std::boolalpha << pairPerm.first << " " << pairPerm.second << std::endl;

		std::vector<int> vecIdx = {0,1,2};
		while(1)
		{
			std::cout << "eps(" << vecIdx[0] << vecIdx[1] << vecIdx[2] << ") = "
				<< tl::epsilon_tensor(vecIdx)
				<< std::endl;
			if(!std::next_permutation(vecIdx.begin(), vecIdx.end()))
				break;
		}

		std::cout << "eps(000) = " << tl::epsilon_tensor({0,0,0}) << std::endl;
		std::cout << "eps(002) = " << tl::epsilon_tensor({0,0,2}) << std::endl;
		std::cout << "eps(110) = " << tl::epsilon_tensor({1,1,0}) << std::endl;
		std::cout << "eps(111) = " << tl::epsilon_tensor({1,1,1}) << std::endl;
		std::cout << "eps(122) = " << tl::epsilon_tensor({1,2,2}) << std::endl;
	}

	return 0;
}
