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

// gcc -o traits traits.cpp -lstdc++ -std=c++14

#include <iostream>
#include <boost/units/systems/si.hpp>
#include "../log/debug.h"
#include "../helper/traits.h"


template<typename T> using t_arr = std::array<T, 3>;


template<typename ...Args> using t_fkt2 = int(*)(Args...);

template<typename t_arg>
using t_fkt_vararg_tst = t_arg(*)(typename std::remove_reference<decltype(((t_arg*)0)[0])>::type);

/*
template<typename t_arg, std::size_t ...idx>
using t_fkt_vararg_impl = t_arg(*)(typename std::remove_reference<decltype(((t_arg*)0)[idx])>::type...);

template<typename t_arg, std::size_t ...idx>
static t_fkt_vararg_impl<t_arg, idx...> tstfkt_vararg(const tl::integer_sequence<std::size_t, idx...>&)
{ return 0; }

template<typename t_arg, std::size_t iNumArgs>
using t_fkt_vararg = decltype(tstfkt_vararg<t_arg>(tl::make_integer_sequence<std::size_t, iNumArgs>()));
*/


int main()
{
	std::cout << short(tl::get_scalar_type<double>::value) << std::endl;
	std::cout << tl::get_typename<tl::get_scalar_type<double>::value_type>() << std::endl;
	std::cout << std::endl;

	std::cout << short(tl::get_scalar_type<boost::units::dimensionless_quantity<boost::units::si::system, double>>::value) << std::endl;
	std::cout << tl::get_typename<tl::get_scalar_type<boost::units::dimensionless_quantity<boost::units::si::system, double>>::value_type>() << std::endl;
	std::cout << std::endl;

	std::vector<int> vec({2,4,6});
	auto fkt = [](int a, int b, int c) -> int{ return a*(b+c); };
	std::cout << tl::get_typename<decltype(fkt)>() << std::endl;
	std::cout << tl::call<3>(fkt, vec) << std::endl;

	t_arr<int> arr({2,4,6});
	std::cout << tl::call<arr.size(), int(int,int,int), int, t_arr>
		([](int a, int b, int c)->int{ return a*(b+c); }, arr) << std::endl;
	std::cout << std::endl;


	using t_fkt = int(*)(int, int);
	std::cout << tl::get_typename<t_fkt>() << std::endl;
	std::cout << tl::get_typename<t_fkt2<int,int,int>>() << std::endl;
	std::cout << tl::get_typename<t_fkt_vararg_tst<int>>() << std::endl;
	//std::cout << tl::get_typename<t_fkt_vararg_impl<int, 1, 2, 3>>() << std::endl;
	std::cout << tl::get_typename<tl::t_fkt_vararg<float, 10>>() << std::endl;

	return 0;
}
