/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
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
