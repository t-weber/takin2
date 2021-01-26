/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o cl_dev cl_dev.cpp -std=c++11 -lOpenCL -lstdc++

#include "../cl/cl.h"
#include <iostream>

int main()
{
	cl::Platform plat;
	cl::Device dev;
	if(!tl::get_best_cl_dev<float>(plat, dev))
	{
		std::cerr << "Cannot get devices." << std::endl;
		return -1;
	}

	std::cout << "Platform: " << plat.getInfo<CL_PLATFORM_NAME>() << std::endl;
	std::cout << "Device: " << dev.getInfo<CL_DEVICE_NAME>() << std::endl;

	return 0;
}
