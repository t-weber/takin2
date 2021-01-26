/**
 * units test
 * @author Tobias Weber <tweber@ill.fr>
 * @date jun-20
 * @license GPLv3, see 'LICENSE' file
 */

#define BOOST_TEST_MODULE La1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/units.h"

using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_equals, t_real, t_types)
{
	std::cout << "kB = " << tl2::kB<t_real> << std::endl;
	std::cout << "muB = " << tl2::muB<t_real> << std::endl;

	std::cout << "kB = "
		<< static_cast<t_real>(tl2::kB<t_real>/tl2::meV<t_real>*tl2::kelvin<t_real>)
		<< std::endl;
	std::cout << "muB = "
		<< static_cast<t_real>(tl2::muB<t_real>/tl2::meV<t_real>*tl2::tesla<t_real>)
		<< std::endl;
}
