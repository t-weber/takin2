/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 10-jun-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -I.. -I/usr/include/libqhullcpp -DUSE_QHULL -o hull1 hull1.cpp -lqhull_r -lqhullcpp
 */

#define BOOST_TEST_MODULE Hull1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>

#include "libs/math20.h"
using namespace tl2_ops;


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_equals, t_real, t_types)
{
	using t_vec = tl2::vec<t_real, std::vector>;

	std::vector<t_vec> vecs =
	{
		{1, 1, -1},
		{1, -1, -1},
		{-1, 1, -1},
		{-1, -1, -1},

		{1, 1, 1},
		{1, -1, 1},
		{-1, 1, 1},
		{-1, -1, 1},
	};


	auto [hull, norms, dists] = tl2_qh::get_convexhull<t_vec>(vecs);

	for(std::size_t faceidx=0; faceidx<hull.size(); ++faceidx)
	{
		const auto& face = hull[faceidx];
		std::cout << "face " << faceidx << ":\n";

		for(const auto& vert : face)
			std::cout << "vertex: " << vert << "\n";
		std::cout << "normal: " << norms[faceidx] << "\n";

		std::cout << std::endl;
	}


	BOOST_TEST(hull.size() == 6*2);
	for(const auto& face : hull)
		BOOST_TEST(face.size() == 3);

	std::cout << "--------------------------------------------------------------------------------" << std::endl;
}
