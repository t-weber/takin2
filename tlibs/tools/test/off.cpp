/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */
// gcc -o off off.cpp -std=c++11 -lstdc++

#include "../file/off.h"

int main(int argc, char** argv)
{
	if(argc < 2)
	{
		std::cerr << "No OFF file given." << std::endl;
		return -1;
	}

	tl::Off3d<double> off;
	if(!off.Load(argv[1]))
	{
		std::cerr << "Error loading " << argv[1] << "." << std::endl;
		return -1;
	}

	off.Optimise(0.01);

	if(!off.Save("/tmp/tst.off"))
	{
		std::cerr << "Error saving OFF." << std::endl;
		return -1;
	}

	return 0;
}
