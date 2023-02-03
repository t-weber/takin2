/**
 * tlibs test file
 * @author Tobias Weber <tweber@ill.fr>
 * @license GPLv2 or GPLv3
 * @date 20-nov-2019
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
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
