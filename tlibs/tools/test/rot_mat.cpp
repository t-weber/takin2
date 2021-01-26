/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
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
