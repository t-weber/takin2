/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I. -o kd3 kd3.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include "../math/kd.h"

using namespace tl;

double fkt(double x)
{
	return std::cos(x*x);
}

int main()
{
	std::list<std::vector<double>> lst;
	for(double dx = 0.; dx<10.; dx+=0.01)
		lst.push_back(std::vector<double>({dx, fkt(dx)}));

	Kd<double> kd(lst);

	while(1)
	{
		std::string str;
		std::cout << "Point: ";
		std::getline(std::cin, str);
		std::istringstream istr(str);
		double dx, dy;
		istr >> dx >> dy;

		const std::vector<double>& vecN = kd.GetNearestNode(std::vector<double>{dx,dy});
		std::cout << "Closest Point: ";
		for(double d : vecN)
			std::cout << d << ", ";
		std::cout << std::endl;
	}

	return 0;
}
