/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I. -o bose bose.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <fstream>
#include "../math/neutrons.h"

int main()
{
	std::ofstream ofstr("bosetst.dat");
	double dT = 80.;

	for(double dE=1.; dE<10.; dE+=0.1)
	{
		double dBose = tl::bose<double>(dE, dT);
		ofstr << dE << "\t" << dBose << "\n";
	}

	return 0;
}
