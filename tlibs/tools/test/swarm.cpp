/**
 * Swarm fitting algorithms
 *
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date Feb-17
 * @license GPLv2 or GPLv3
 *
 * gcc -o swarm swarm.cpp ../math/rand.cpp ../log/log.cpp -lstdc++ -lm -lpthread
 */

#include "../fit/swarm.h"
namespace ublas = tl::ublas;

using t_real = double;
template<class T> using t_vec = ublas::vector<T>;

int main()
{
	t_real x[] = {1., 2., 3., 4., 5.};
	t_real y[] = {2., 4., 6., 8., 10.};
	t_real yerr[] = {0.5, 0.5, 0.5, 0.5, 0.5};

	auto func = [](t_real x, t_real m, t_real t) -> t_real { return m*x + t; };
	tl::FitterLamFuncModel<t_real, 3, decltype(func)> mod(func);

	tl::Unkindness<t_real, t_vec> unk;
	unk.SetData(5, x, y, yerr);
	unk.SetMaxCalls(20);
	unk.Init(100, &mod,
		tl::make_vec({0.,0.}), tl::make_vec({1.,1.}));
	unk.Run();

    return 0;
}
