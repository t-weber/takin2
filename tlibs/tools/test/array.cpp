/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o array array.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include "../helper/array.h"

int main()
{
	int is[] = {1,2,3,4,5};
	tl::wrapper_array<int> a(is, 5);

	for(int i : a)
		std::cout << i << std::endl;

	return 0;
}
