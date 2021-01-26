/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o proc proc.cpp -lstdc++ -lboost_iostreams

#include "../helper/proc.h"
#include <iostream>

int main(int argc, char **argv)
{
	if(argc < 2) return -1;

	tl::PipeProc<> proc(argv[1], 0);
	if(proc.IsReady())
	{
		std::istream& istr = proc.GetIstr();
		while(!istr.eof())
		{
			std::string strLine;
			std::getline(istr, strLine);
			std::cout << strLine << std::endl;
		}
	}

	return 0;
}
