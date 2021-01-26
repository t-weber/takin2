/**
 * tlibs test file
 * @author Tobias Weber <tweber@ill.fr>
 * @license GPLv2 or GPLv3
 * @date 20-nov-2019
 */

// g++ -I../.. -o quat2 quat2.cpp -std=c++11

#include "math/quat.h"
#include "math/linalg.h"
#include <iostream>

using namespace tl;
using T = double;

int main()
{
	ublas::matrix<T> matX = rotation_matrix_3d_x(45./180.*M_PI);
	ublas::matrix<T> matY = rotation_matrix_3d_y(66./180.*M_PI);

	math::quaternion<T> quatX = rot3_to_quat(matX);
	math::quaternion<T> quatY = rot3_to_quat(matY);

	math::quaternion<T> quat = slerp(quatX, quatY, 0.25);

	ublas::matrix<T> mat = quat_to_rot3(quat);


	std::cout << "matX = " << matX << std::endl;
	std::cout << "matY = " << matY << std::endl;
	std::cout << "mat  = " << mat << std::endl;

	return 0;
}
