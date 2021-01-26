/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o scatangles test/scatangles.cpp -lm -lstdc++ -std=c++11

#include <iostream>
#include "../math/neutrons.h"

int main()
{
	using T = float;

	T dKi = 1.4;
	T dKf = 1.4;
	T dQ = 2.5;

	while(1)
	{
		std::cout << "ki = "; std::cin >> dKi;
		std::cout << "kf = "; std::cin >> dKf;
		std::cout << "Q = "; std::cin >> dQ;

		tl::t_wavenumber_si<T> ki = dKi / tl::get_one_angstrom<T>();
		tl::t_wavenumber_si<T> kf = dKf / tl::get_one_angstrom<T>();
		tl::t_wavenumber_si<T> Q = dQ / tl::get_one_angstrom<T>();

		T kiQ = tl::get_angle_ki_Q(ki, kf, Q, 1, 0) / tl::get_one_radian<T>();
		T kfQ = tl::get_angle_kf_Q(ki, kf, Q, 1, 0) / tl::get_one_radian<T>();

		std::cout << "kiQ = " << tl::r2d(kiQ) << std::endl;
		std::cout << "kfQ = " << tl::r2d(kfQ) << std::endl;

		std::cout << std::endl;
	}

	return 0;
}
