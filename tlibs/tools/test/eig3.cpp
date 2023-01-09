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

// gcc -I/usr/include/lapacke -o eig3 eig3.cpp ../math/linalg2.cpp ../log/log.cpp -lstdc++ -lm -llapacke -llapack -std=c++11

#define TLIBS_INC_HDR_IMPLS
#include "../math/linalg.h"
#include "../math/linalg2.h"

namespace ublas = boost::numeric::ublas;
using T = double;

int main()
{
	ublas::matrix<T> M = tl::make_mat({
		{1.1,-2.2,3.3},
		{-2.2,5.5,8.8},
		{3.3,8.8,-9.9}	});

	ublas::matrix<T> M_org = M;
	std::cout << M << std::endl;

	std::vector<ublas::vector<T>> evecs;
	std::vector<T> evals;
	tl::eigenvec_sym<T>(M, evecs, evals);
	for(int i=0; i<evals.size(); ++i)
		std::cout << "eval: " << evals[i] <<
		", evec: " << (evecs[i]/ublas::norm_2(evecs[i])) <<
		", len: " << ublas::norm_2(evecs[i]) << std::endl;
	std::cout << std::endl;

	double eval;
	ublas::vector<double> evec;
	tl::eigenvec_dominant_sym(M, evec, eval);
	std::cout << "dominant evec: " << evec << std::endl;
	std::cout << "dominant eval: " << eval << std::endl;
	tl::eigenvec_least_dominant_sym(M, evec, eval);
	std::cout << "least dominant evec: " << evec << std::endl;
	std::cout << "least dominant eval: " << eval << std::endl;
	std::cout << std::endl;

	return 0;
}
