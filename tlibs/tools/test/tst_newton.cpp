/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

#include "../math/numint.h"
#include <iostream>

int main()
{
	double dRes = tl::newton<double>(
		[](double x)->double { return x*x*x*x + x*x*x + x*x + x + 1.; },
		[](double x)->double { return 4*x*x*x + 3*x*x + 2*x + 1.; },
		0., 512, 0.0001);
	std::cout << dRes << std::endl;

	return 0;
}
