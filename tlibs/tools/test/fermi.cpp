/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I. -o fermi fermi.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <fstream>
#include "../math/neutrons.h"

using T = double;


int main()
{
	std::ofstream ofstr("fermitst.dat");
	T dMu = 5.;

	for(T dT=0.; dT<20.; dT+=5.)
	{
		for(T dE=0.; dE<10.; dE+=0.1)
		{
			T n = tl::fermi<T>(dE, dMu, dT);
			ofstr << dE << "\t" << n << "\n";
		}
		ofstr << "e\n";
	}

	return 0;
}
