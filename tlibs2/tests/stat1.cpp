/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 20-aug-20
 * @license GPLv3, see 'LICENSE' file
 */

#define BOOST_TEST_MODULE Stat1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>
#include <boost/math/quaternion.hpp>

#include "libs/math20.h"
using namespace tl2_ops;


using t_types = std::tuple<long double, double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_equals, t_real, t_types)
{
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;

	std::vector<t_vec> vecs = {{
		tl2::create<t_vec>({-1, 60}),
		tl2::create<t_vec>({20, 5}),
		tl2::create<t_vec>({-3, 40}),
		tl2::create<t_vec>({40, 3}),
		tl2::create<t_vec>({-5, 20}),
		tl2::create<t_vec>({60, 1}),
	}};

	auto [cov, cor] = tl2::covariance<t_mat, t_vec>(vecs);
	std::cout << "cov = " << cov << std::endl;
	std::cout << "cor = " << cor << std::endl;

	//BOOST_TEST(tl2::equals(vec, vec2, 1e-5));
}
