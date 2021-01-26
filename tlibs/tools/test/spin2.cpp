/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 * clang -I/usr/include/lapacke -o spin2 ../test/spin2.cpp ../log/log.cpp -lstdc++ -std=c++11 -lm
 */

#include "../phys/spin.h"
#include <iostream>

using t_real = double;
using namespace tl;

int main()
{
	t_real J, j1, j2, m1, m2;
	std::cout << "J = "; std::cin >> J;
	std::cout << "j1 = "; std::cin >> j1;
	std::cout << "j2 = "; std::cin >> j2;
	std::cout << "m1 = "; std::cin >> m1;
	std::cout << "m2 = "; std::cin >> m2;

	std::cout << CG_coeff(J, j1, j2, m1, m2) << std::endl;
	return 0;
}
