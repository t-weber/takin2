/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

#include <iostream>

#define HAS_COMPLEX_ERF
#include "../math/math.h"

int main()
{
	std::cout << tl::voigt_model<double>(0., 0., 1., 1., 1., 0.) << std::endl;
	return 0;
}
