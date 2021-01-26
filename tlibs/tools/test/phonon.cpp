/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I. -o phonon test/phonon.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <fstream>
#include "../math/neutrons.h"


int main()
{
	std::ofstream ofstr("phonontst.gpl");

	double dT = 120.;
	double dE0 = 2.;
	double dHWHM = 1.;
	double dAmp = 1.;
	double dOffs = 0.;

	ofstr << "plot \"-\" pointtype 7\n";
	for(double dE=-5.; dE<5.; dE+=0.1)
	{
		double dS = tl::DHO_model<double>(dE, dT, dE0, dHWHM, dAmp, dOffs);
		ofstr << dE << "\t" << dS << "\n";
	}
	ofstr << "end\n";

	return 0;
}
