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

// gcc -I/usr/include/lapacke -o eig eig.cpp ../math/linalg2.cpp ../log/log.cpp -lstdc++ -lm -llapacke -llapack -std=c++11

#include "../math/linalg.h"
#include "../math/linalg2.h"

namespace ublas = boost::numeric::ublas;

int main()
{
	ublas::matrix<double> M = tl::make_mat({{1.1,-2.2,3.3},
						{-2.2,5.5,8.8},
						{3.3,8.8,-9.9}});
	ublas::matrix<double> M_org = M;
	std::cout << M << std::endl;

	std::vector<ublas::vector<double>> evecs;
	std::vector<double> evals;
	tl::eigenvec_sym<double>(M, evecs, evals);
	for(int i=0; i<evals.size(); ++i)
		std::cout << "eval: " << evals[i] << ", evec: " << (evecs[i]/ublas::norm_2(evecs[i])) << std::endl;
	std::cout << std::endl;

	ublas::matrix<double> I = ublas::identity_matrix<double>(3);

	for(int i=0; i<10; ++i)
	{
		ublas::matrix<double> Q, R;
		if(tl::qr_decomp(M, Q, R))
		{
			std::cout << "Q = " << Q << std::endl;
			std::cout << "R = " << R << std::endl;

			//Q = tl::norm_col_vecs(Q);
			M = ublas::prod(R, Q);
			I = ublas::prod(I, Q);

			/*ublas::matrix<double> QT = ublas::trans(Q);
			QT = ublas::prod(QT, M_org);
			QT = ublas::prod(QT, Q);
			std::cout << "Q^t M Q = " << QT << std::endl;*/

			std::cout << "M = " << M << std::endl;
			std::cout << "I = " << I << std::endl;
			std::cout << std::endl;
		}
		else
		{
			std::cerr << "Cannot determine QR." << std::endl;
			break;
		}
	}
	return 0;
}
