/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

#include "../file/file.h"

int main()
{
	std::vector<std::string> vecPaths =
	{
		"/usr/include/", 
		"/usr/local/include/", 
	};


	bool bFound;
	std::string strFile;

	std::tie(bFound, strFile) = tl::find_file(vecPaths, std::string("stdio.h"));
	std::cout << "found: " << std::boolalpha << bFound << ", ";
	std::cout << strFile << std::endl;

	std::tie(bFound, strFile) = tl::find_file(vecPaths, std::string("./mpc.h"));
	std::cout << "found: " << std::boolalpha << bFound << ", ";
	std::cout << strFile << std::endl;

	std::tie(bFound, strFile) = tl::find_file(vecPaths, std::string("/mpc.h"));
	std::cout << "found: " << std::boolalpha << bFound << ", ";
	std::cout << strFile << std::endl;

	std::tie(bFound, strFile) = tl::find_file(vecPaths, std::string("nosuchfile"));
	std::cout << "found: " << std::boolalpha << bFound << ", ";
	std::cout << strFile << std::endl;

	std::tie(bFound, strFile) = tl::find_file(vecPaths, std::string("/opt/X11/include/png.h"));
	std::cout << "found: " << std::boolalpha << bFound << ", ";
	std::cout << strFile << std::endl;

	return 0;
}
