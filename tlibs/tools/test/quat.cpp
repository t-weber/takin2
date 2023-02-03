/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
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

// gcc -o quat quat.cpp -std=c++11 -lstdc++ -lm

#include "../math/quat.h"
#include <iostream>

using namespace tl;
typedef float T;

int main()
{
	T angle = 1.23;
	ublas::matrix<T> mat = rotation_matrix_3d_x(angle);
	std::cout << "mat = " << mat << std::endl;

	math::quaternion<T> quat = rot3_to_quat(mat);
	std::cout << "quat = " << quat << std::endl;

	mat = quat_to_rot3(quat);
	std::cout << "mat = " << mat << std::endl;

	std::cout << quat_to_cmat(quat) << std::endl;


	ublas::vector<T> vec = make_vec({1., 2., 3.});
	std::cout << ublas::prod(mat, vec) << std::endl;
	std::cout << quat_vec_prod(quat, vec) << std::endl;



	ublas::matrix<T> mat2 = rotation_matrix_3d_y(angle);
	math::quaternion<T> quat2 = rot3_to_quat(mat2);

	std::cout << "q*q2 = " << quat*quat2 << std::endl;
	auto vecEuler = quat_to_euler(quat*quat2);
	std::cout << "euler xyz: ";
	for(T ang : vecEuler)
		std::cout << ang/M_PI*180. << " ";
	std::cout << std::endl;
	std::cout << "q from euler xyz: " 
		<< euler_to_quat(vecEuler[0],vecEuler[1],vecEuler[2]) 
		<< std::endl;
	/*auto vecEuler2 = quat_to_euler_zxz(quat*quat2);
	std::cout << "euler zxz: ";
	for(T ang : vecEuler2)
		std::cout << ang/M_PI*180. << " ";
	std::cout << std::endl;
	std::cout << "q from euler zxz: "
		<< euler_to_quat_zxz(vecEuler2[0],vecEuler2[1],vecEuler2[2]) 
		<< std::endl;*/


	std::cout << "prod1: " << quat_prod(quat, quat2) << std::endl;
	std::cout << "prod2: " << quat * quat2 << std::endl;
	return 0;
}
