/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
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
