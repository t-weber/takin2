/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// clang -o varg varg.cpp -std=c++14 -lstdc++ -lm
// ideas for future tl::ThreadPool class with args


#include <vector>
#include <functional>
#include <tuple>
#include <utility>
#include <iostream>
#include "../helper/misc.h"
#include "../log/debug.h"


template<class ...t_args>
using t_vecArgs = std::vector<std::tuple<t_args...>>;

template<class t_ret, class ...t_args>
using t_vecFuncs = std::vector<std::function<t_ret(t_args...)>>;


template<class t_ret, class ...t_args>
void fkttst(t_ret (&&func)(t_args...), t_args&&... args)
{
	t_vecFuncs<t_ret, t_args...> vecFuncs = {func};
	t_vecArgs<t_args...> vecArgs = {std::tuple<t_args...>(args...)};
	constexpr int iNumArgs = std::tuple_size<typename decltype(vecArgs)::value_type>::value;

	std::cout << tl::get_typename<decltype(vecFuncs)>() << std::endl;
	std::cout << iNumArgs << " args: " << tl::get_typename<decltype(vecArgs)>() << std::endl;
	std::cout << std::endl;


	for(std::size_t i=0; i<vecFuncs.size(); ++i)
	{
		typename decltype(vecFuncs)::value_type& func = vecFuncs[i];
		typename decltype(vecArgs)::value_type& arg = vecArgs[i];

		auto idx = std::index_sequence_for<t_args...>();

		auto callfkt = [&func, &idx, &arg](auto&& idx) -> t_ret
		{ return func(std::get<t_args>(arg)...); };

		t_ret ret = callfkt(idx);

		std::cout << "Func 1: " << std::endl;
		std::cout << tl::get_typename<decltype(func)>() << std::endl;
		std::cout << tl::get_typename<decltype(arg)>() << std::endl;
		std::cout << "Index seq: " << tl::get_typename<decltype(idx)>() << std::endl;
		std::cout << "Return value: " << ret << std::endl;
		std::cout << std::endl;
	}
}


int tst0(int i, int j) { return i*j; }

int main()
{
	fkttst<int, int, int>(tst0, 2, 3);
}
