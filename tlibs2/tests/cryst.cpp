/**
 * crystal matrix text
 * @author Tobias Weber <tweber@ill.fr>
 * @date feb-19
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -o cryst cryst.cpp
 * g++ -std=c++20 -DUSE_LAPACK -I/usr/include/lapacke -I/usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -o cryst cryst.cpp -llapacke
 */

#define BOOST_TEST_MODULE Xtal Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;


#include <iostream>
#include <vector>

#include "libs/math20.h"
using namespace tl2_ops;


#ifdef USE_LAPACK
	using t_types = std::tuple<double, float>;
#else
	using t_types = std::tuple<long double, double, float>;
#endif
BOOST_AUTO_TEST_CASE_TEMPLATE(test_xtal, t_real, t_types)
{
	//using t_cplx = std::complex<t_real>;
	//using t_vec = std::vector<t_real>;
	using t_mat = tl2::mat<t_real, std::vector>;
	//using t_vec_cplx = std::vector<t_cplx>;
	//using t_mat_cplx = tl2::mat<t_cplx, std::vector>;

	auto A = tl2::A_matrix<t_mat, t_real>(3., 4., 5., 80./180.*M_PI, 100./180.*M_PI, 60./180.*tl2::pi<t_real>);
	auto B = tl2::B_matrix<t_mat, t_real>(3., 4., 5., 80./180.*M_PI, 100./180.*M_PI, 60./180.*tl2::pi<t_real>);
	auto [B2, ok] = tl2::inv<t_mat>(A);
	B2 = 2.*tl2::pi<t_real> * tl2::trans<t_mat>(B2);

	std::cout << "A  = " << A << std::endl;
	std::cout << "B  = " << B << std::endl;
	std::cout << "B2 = " << B2 << std::endl;

	BOOST_TEST(ok);
	BOOST_TEST(tl2::equals(B, B2, std::numeric_limits<t_real>::epsilon()*1e2));
	std::cout << std::endl;
}
