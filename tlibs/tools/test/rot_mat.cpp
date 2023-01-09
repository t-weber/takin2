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

// gcc -o rot_mat rot_mat.cpp -std=c++11 -lstdc++ -lm

#include "../math/quat.h"
#include <iostream>

using namespace tl;
using T = double;
using t_vec = ublas::vector<T>;
using t_mat = ublas::matrix<T>;


int main()
{
	t_vec vec0 = make_vec<t_vec>({0,0,1});
	t_vec vec1 = make_vec<t_vec>({1,1,1});

	T vx, vy, vz;
	std::cout << "vec0 = ";
	std::cin >> vx >> vy >> vz;
	vec0[0] = vx; vec0[1] = vy; vec0[2] = vz;

	std::cout << "vec1 = ";
	std::cin >> vx >> vy >> vz;
	vec1[0] = vx; vec1[1] = vy; vec1[2] = vz;

	std::cout << "Rotating " << vec0 << " into " << vec1 << "." << std::endl;

	auto quat = rotation_quat(vec0, vec1);
	std::cout << "quat = " << quat << std::endl;

	auto mat = quat_to_rot3(quat);
	std::cout << "mat = " << mat << std::endl;
	std::cout << std::endl;

	//std::cout << quat_vec_prod(quat, vec0) << std::endl;
	//std::cout << ublas::prod(mat, vec0) << std::endl;


	T dScale = 1.;
	std::cout << "scale = ";
	std::cin >> dScale;

	t_vec vecBase1 = make_vec<t_vec>({1,0,0});
	t_vec vecBase2 = make_vec<t_vec>({0,1,0});
	t_vec vecBase3 = make_vec<t_vec>({0,0,1});
	std::cout << "quat * [100] * scale = "
		<< quat_vec_prod<t_vec>(quat, vecBase1)*dScale << std::endl;
	std::cout << "quat * [010] * scale = "
		<< quat_vec_prod<t_vec>(quat, vecBase2)*dScale << std::endl;
	std::cout << "quat * [001] * scale = "
		<< quat_vec_prod<t_vec>(quat, vecBase3)*dScale << std::endl;
	std::cout << std::endl;

	return 0;
}
