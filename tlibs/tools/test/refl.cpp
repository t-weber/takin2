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

// gcc -o refl refl.cpp -lstdc++ -std=c++11 -lm
#include "../math/linalg.h"

namespace ublas = boost::numeric::ublas;
using t_real = float;

int main()
{
	ublas::vector<t_real> vec = tl::make_vec({1.,2.,3.});
	ublas::vector<t_real> vecN = tl::make_vec({1.,2.,0.});
	ublas::vector<t_real> vecRefl = tl::reflection(vec, vecN);

	std::cout << "vec: " << vec << std::endl;
	std::cout << "norm: " << vecN << std::endl;
	std::cout << "refl: " << vecRefl << std::endl;


	ublas::matrix<t_real> M = tl::make_mat({{1.1,-2.2,3.3},{-4.4,5.5,6.6}, {7.7,8.8,-9.9}});
	std::cout << M << std::endl;
	std::cout << tl::insert_unity(M, 2) << std::endl;


	std::cout << "--------------------------------------" << std::endl;
	std::cout << "qr via householder algo:" << std::endl;
	ublas::matrix<t_real> Q, R;
	if(tl::qr_decomp(M, Q, R))
	{
		std::cout << "M = " << M << std::endl;
		std::cout << "Q = " << Q << std::endl;
		std::cout << "R = " << R << std::endl;
		std::cout << "QR = " << ublas::prod(Q,R) << std::endl;
	}

	std::cout << "--------------------------------------" << std::endl;
	std::cout << "qr via gram-schmidt algo:" << std::endl;
	ublas::matrix<t_real> Q0, R0;
	if(tl::qr_decomp_gs(M, Q0, R0))
	{
		std::cout << "M = " << M << std::endl;
		std::cout << "Q = " << Q0 << std::endl;
		std::cout << "R = " << R0 << std::endl;
		std::cout << "QR = " << ublas::prod(Q0,R0) << std::endl;
	}

	return 0;
}
