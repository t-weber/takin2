/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o ostream ostream.cpp ../log/log.cpp -lstdc++ -std=c++11 -lboost_iostreams

#include "../file/comp.h"
#include <memory>

int main()
{
	std::ofstream ofstr("tst.txt.gz");
	std::shared_ptr<std::ostream> pOstr = tl::create_comp_ostream(ofstr, tl::Compressor::GZ);

	for(int i=0; i<1000; ++i)
		(*pOstr) << i << "\n";

	return 0;
}
