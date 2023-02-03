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

// gcc -I. -o mat1 mat1.cpp -std=c++11 -lstdc++ -lm


#include "../math/linalg.h"
#include <iostream>

using T = double;
using t_mat = tl::ublas::matrix<T>;
using t_vec = tl::ublas::vector<T>;


int main()
{
	std::vector<t_vec> vecs =
	{
		tl::make_vec({1., 2., 3.}),
		tl::make_vec({2., 1., 4.}),
		tl::make_vec({5., 1., 3.}),
		tl::make_vec({2., 2., 3.}),
		tl::make_vec({6., 4., 2.}),
	};

	auto vecMean = tl::mean_value(vecs);
	for(t_vec& vec : vecs)
		vec -= vecMean;

	t_mat mat = tl::column_matrix(vecs);

	t_mat matCov1, matCorr1;
	std::tie(matCov1, matCorr1) = tl::covariance(vecs);
	t_mat matCov2 = tl::ublas::prod(mat, tl::ublas::trans(mat));

	std::cout << "cov(...) = " << matCov1 << std::endl;
	std::cout << "corr(...) = " << matCorr1 << std::endl;
	std::cout << "AA^t =     " <<  matCov2/matCov2(0,0)*matCov1(0,0) << std::endl;

	return 0;
}
