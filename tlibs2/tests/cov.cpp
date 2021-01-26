/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date nov-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -I.. -o cov cov.cpp
 */

#define BOOST_TEST_MODULE Cov1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/math20.h"
using namespace tl2_ops;


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_equals, t_real, t_types)
{
	using t_vec = std::vector<t_real>;
	using t_mat = tl2::mat<t_real, std::vector>;

	std::vector<t_vec> vecs{{
		tl2::create<t_vec>({1, 2, 3}),
		tl2::create<t_vec>({3, 2, 1}),
		tl2::create<t_vec>({5, -2, 0.5}),
	}};

	std::vector<t_real> probs{{ 1., 1., 0.5 }};

	auto [cov, corr] = tl2::covariance<t_mat, t_vec>(vecs, &probs);
	std::cout << "cov  = " << cov << std::endl;
	std::cout << "corr = " << corr << std::endl;

	// comparing with: np.cov([[1,2,3], [3,2,1], [5,-2,0.5]], rowvar=False, aweights=[1,1,0.5], ddof=0)
	const t_mat cov_cmp = tl2::create<t_mat>({
		 2.24, -1.92, -1.52,
		-1.92,  2.56,  0.96,
		-1.52,  0.96,  1.16
	});
	BOOST_TEST(tl2::equals(cov, cov_cmp, 1e-5));
}
