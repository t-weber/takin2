/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I. -o mat0 mat0.cpp -std=c++11 -lstdc++ -lm

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
