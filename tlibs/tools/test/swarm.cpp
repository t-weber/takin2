/**
 * Swarm fitting algorithms
 *
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date Feb-17
 * @license GPLv2 or GPLv3
 *
 * gcc -o swarm swarm.cpp ../math/rand.cpp ../log/log.cpp -lstdc++ -lm -lpthread
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
