/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o mag mag.cpp ../math/rand.cpp ../log/log.cpp -lstdc++ -lm

#include <fstream>
#include "../math/mag.h"

using T = double;

int main()
{
	T J = 1;

	std::ofstream ofstr("brill.dat");

	for(T x=0; x<5.; x+=0.01)
	{
		// x ~ B/T
		ofstr << x << "\t" << (1./x) << "\t" << tl::brillouin<T>(J, x);
		ofstr << "\n";
	}

	return 0;
}
