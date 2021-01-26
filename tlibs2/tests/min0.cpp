/**
 * minimisation test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 27-aug-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++-10 -std=c++20 -I.. -o min min0.cpp ../libs/log.cpp -lMinuit2
 */

#include <string>
#include <iostream>

#define BOOST_TEST_MODULE Min Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include "libs/fit.h"


using t_types_real = std::tuple<double, float>;


BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_real, t_real, t_types_real)
{
	std::string func{"(x-5.67)^2"};

	std::vector<std::string> params{{"x"}};
	std::vector<t_real> vals{{0.}};
	std::vector<t_real> errs{{0.1}};
	std::vector<bool> fixed{false};

	bool ok = tl2::minimise_expr(func, params, vals, errs, &fixed);

	BOOST_TEST(ok);
	BOOST_TEST(tl2::equals<t_real>(vals[0], 5.67, 1e-3));
}
