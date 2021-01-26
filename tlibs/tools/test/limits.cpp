/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I. -o limits limits.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <limits>

template<typename T>
void print_limits()
{
	std::cout << "epsilon: " << std::numeric_limits<T>::epsilon() << std::endl;
	std::cout << "lowest: " << std::numeric_limits<T>::lowest() << std::endl;

	std::cout << "min: " << std::numeric_limits<T>::min() << std::endl;
	std::cout << "max: " << std::numeric_limits<T>::max() << std::endl;

	std::cout << "digits: " << std::numeric_limits<T>::digits << std::endl;
	std::cout << "digits10: " << std::numeric_limits<T>::digits10 << std::endl;
	std::cout << "max_digits10: " << std::numeric_limits<T>::max_digits10 << std::endl;
}


int main()
{
	std::cout << "Limits for double:" << std::endl;
	print_limits<double>();
	std::cout << std::endl;

	std::cout << "Limits for float:" << std::endl;
	print_limits<float>();
	return 0;
}
