/**
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-17
 * @license GPLv2 or GPLv3
 */

#include "../phys/neutrons.h"


int main()
{
	auto A = tl::get_one_angstrom<double>();
	auto m = tl::get_one_meter<double>();
	auto deg = tl::get_one_deg<double>();

	std::cout << tl::vsel_freq(2.36*A, 0.1*m, 50.*deg) << std::endl;

	return 0;
}
