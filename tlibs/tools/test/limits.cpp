/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
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
