/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o fs fs.cpp -lstdc++ -lboost_system -lboost_filesystem -std=c++11

#include <iostream>
#include "../file/file.h"

int main()
{
	std::cout << tl::get_file_size(std::string("test/fs.cpp")) << std::endl;

	std::vector<std::string> vec = tl::get_all_files<0>("/home/tweber/tmp");
	std::vector<std::string> vecRec = tl::get_all_files<1>("/home/tweber/tmp");

	std::cout << "\nnon-recursive: " << std::endl;
	for(const auto& file : vec)
		std::cout << file << std::endl;
	std::cout << "\n\nrecursive: " << std::endl;
	for(const auto& file : vecRec)
		std::cout << file << std::endl;
	return 0;
}
