/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -DNO_LAPACK -o geo geo.cpp ../log/log.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include "../math/geo.h"

using T = double;
using t_vec = typename tl::Plane<T>::t_vec;

int main()
{
	tl::Plane<T> plane1(
		tl::make_vec({1., 1., 0.}),	// vec0
		tl::make_vec({0., 1.5, 1.}),	// dir0
		tl::make_vec({1., 0., 1.})	// dir1
		);

	tl::Plane<T> plane2(
		tl::make_vec({0., 2., 0.}),	// vec0
		tl::make_vec({1., 0., 0.}),	// dir0
		tl::make_vec({0., 0., -1.})	// dir1
		);

	tl::Plane<T> plane3(
		tl::make_vec({1., 1., 3.}),	// vec0
		tl::make_vec({1., 1., 0.}),	// dir0
		tl::make_vec({0., 2., 0.5})	// dir1
		);

	t_vec pt;
	plane1.intersect(plane2, plane3, pt);
	std::cout << pt << std::endl;


	// alternative method
	tl::Line<T> line;
	plane1.intersect(plane2, line);

	T t;
	line.intersect(plane3, t);
	std::cout << line(t) << std::endl;
}
