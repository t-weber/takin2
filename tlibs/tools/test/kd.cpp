/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I. -o kd kd.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <vector>
#include <list>
#include "../math/kd.h"

int main()
{
	std::list<std::vector<int>> lst =
	{
		{1,9},
		{2,8},
		{3,7},
		{4,6},
		{5,5},
		{6,4},
		{7,3},
		{8,2},
		{9,1},
		{10,0},
		{20,15}
	};

	tl::Kd<int> kd(lst);
	kd.GetRootNode()->print(std::cout);
	std::cout << std::endl;

	const std::vector<int>& vecN = kd.GetNearestNode(std::vector<int>{20,20});
	std::cout << "nearest: ";
	for(int i : vecN)
		std::cout << i << ", ";
	std::cout << std::endl;
}
