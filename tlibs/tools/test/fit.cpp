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

// gcc -std=c++14 -I/usr/include/root -o fit fit.cpp ../log/log.cpp -lstdc++ -L/usr/lib64/root -lMinuit2

#include "../fit/minuit.h"
#include "../log/debug.h"
#include <boost/type_traits/function_traits.hpp>

using T = tl::t_real_min;

// test instantiation
template class tl::Chi2Function_nd<T, std::vector>;


template<class t_func, class... t_args>
void tst(t_func&& func)
{
	// these do not work with lambdas...
	std::cout << sizeof...(t_args) << std::endl;
	std::cout << boost::function_traits<t_func(t_args...)>::arity << std::endl;
}

int main()
{
	tst([](T x, T m)->T { return m*x; });

	std::vector<T> vecX{0., 1., 2., 3., 4.};
	std::vector<T> vecY{0., 1.2, 1.8, 3.1, 4.2};
	std::vector<T> vecYErr{0.5, 0.5, 0.5, 0.5, 0.5};

	std::vector<std::string> vecParamNames = { "m", "t" };
	std::vector<T> vecVals = { 0.5, 0. };
	std::vector<T> vecErrs = { 0.5, 0. };

	bool bOk = tl::fit<3>([](T x, T m, T t)->T { return m*x + t; },
		vecX, vecY, vecYErr,
		vecParamNames, vecVals, vecErrs);

	std::cout << "Fit: valid = " << bOk 
		<< ", m = " << vecVals[0] << " +- " << vecErrs[0]
		<< ", t = " << vecVals[1] << " +- " << vecErrs[1]
		<< std::endl;

	return 0;
}
