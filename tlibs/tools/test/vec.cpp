/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o vec vec.cpp ../log/log.cpp -lstdc++ -lm -std=c++11

#define TLIBS_USE_GLOBAL_OPS
#include "../math/linalg.h"
#include "../math/linalg_ops.h"

namespace ublas = boost::numeric::ublas;

/*template<class t_mat=ublas::matrix<double>>
t_mat operator*(const t_mat& mat0, const t_mat& mat1)
{
	return ublas::prod(mat0, mat1);
}

template<class t_mat=ublas::matrix<double>,
	class t_vec=ublas::vector<double>>
t_vec operator*(const t_mat& mat, const t_vec& vec)
{
	return ublas::prod(mat, vec);
}

template<class t_vec=ublas::vector<double>>
typename t_vec::value_type operator*(const t_vec& vec0, const t_vec& vec1)
{
	return ublas::inner_prod(vec0, vec1);
}*/

int main()
{
	ublas::matrix<double> M0 = tl::make_mat({{1.1,-2.2,3.3},
						{-2.2,5.5,8.8},
						{3.3,8.8,-9.9}});

	ublas::matrix<double> M1 = tl::make_mat({{1.1,-2.2,3.3},
						{-2.2,5.5,8.8},
						{3.3,8.8,-9.9}});

	ublas::vector<double> v = tl::make_vec({1.,2.,3.});

	std::cout << M0*M1 << std::endl;
	std::cout << M0*v << std::endl;
	std::cout << v*v << std::endl;

	return 0;
}
