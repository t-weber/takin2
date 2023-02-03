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

// g++ -I.. -o mat0 mat0.cpp -std=c++11

#include "../math/linalg.h"
#include <iostream>

using namespace tl;

int main()
{
	//std::cout << std::numeric_limits<double>::min() << std::endl;
	//std::cout << -std::numeric_limits<double>::max() << std::endl;

	ublas::matrix<double> mat1 = make_mat({{1., 2.}, {2.5, 3.}});
	ublas::matrix<double> mat2 = make_mat({{1., 2.}, {2., 3.}});
	ublas::matrix<double> matA = make_mat({{1., 2.}, {3., 4.}});

	std::cout << tl::prod_mm(mat1, mat2) << std::endl;

	std::cout << "row: " << get_row(mat1,0) << std::endl;
	std::cout << "col: " << get_column(mat1,0) << std::endl;
	std::cout << "collen: " << veclen(get_column(mat1,0)) << std::endl;
	std::cout << "inner: " << inner(get_column(mat1,0), get_column(mat1,1)) << std::endl;
	std::cout << "outer: " << outer(get_column(mat1,0), get_column(mat1,1)) << std::endl;
	std::cout << "block: " << block_matrix(mat1, mat2) << std::endl;

	double(*sqrt)(double) = std::sqrt;
	ublas::matrix<double> mat3 = apply_fkt(mat2, std::function<double(double)>(sqrt));

	std::cout << "prod m*v: " << tl::prod_mv(matA, get_column(mat1,1)) << std::endl;
	std::cout << "prod m*v: " << ublas::prod(matA, get_column(mat1,1)) << std::endl;
	std::cout << "prod v*m: " << tl::prod_vm(get_column(mat1,1), matA) << std::endl;
	std::cout << "prod v*m: " << ublas::prod(get_column(mat1,1), matA) << std::endl;

	std::cout << mat3 << std::endl;

	std::cout << is_symmetric(mat1) << std::endl;
	std::cout << is_symmetric(mat2) << std::endl;

	std::cout << get_minmax(mat1).first << std::endl;
	std::cout << get_minmax(mat1).second << std::endl;
	return 0;
}
