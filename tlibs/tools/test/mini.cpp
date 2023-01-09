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

// gcc -std=c++14 -I/usr/include/root -o mini mini.cpp ../log/log.cpp -lstdc++ -L/usr/lib64/root -lMinuit2

#include "../fit/minuit.h"
#include "../log/debug.h"
#include <boost/type_traits/function_traits.hpp>

using T = tl::t_real_min;

int main()
{
	std::vector<std::string> vecParamNames = { "x" };
	std::vector<T> vecVals = { 0.5, };
	std::vector<T> vecErrs = { 0.5, };

	auto fkt = [](T x)->T { return x*x + 10.*x + 10.; };
	bool bOk = tl::minimise<1>(fkt, vecParamNames, vecVals, vecErrs);

	std::cout << "Minumum: valid = " << bOk
		<< ", x = " << vecVals[0] << " +- " << vecErrs[0]
		<< ", f(x) = " << fkt(vecVals[0])
		<< std::endl;

	return 0;
}
