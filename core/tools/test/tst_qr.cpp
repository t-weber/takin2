/**
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

// gcc -I../.. -DUSE_LAPACK -o tst_qr tst_qr.cpp ../../tlibs/math/linalg2.cpp ../../tlibs/log/log.cpp -std=c++11 -lstdc++ -lm -I/usr/include/lapacke -llapacke -llapack

#include "tlibs/math/linalg.h"
#include "tlibs/math/linalg2.h"

using namespace tl;

typedef ublas::vector<double> vec;
typedef ublas::matrix<double> mat;

int main()
{
	mat M(2,4);
	M(0,0)=1; M(0,1)=8; M(0,2)=2; M(0,3)=0;
	M(1,0)=2; M(1,1)=0; M(1,2)=2; M(1,3)=1;
	//M = ublas::trans(M);
	std::cout << "M = " << M << std::endl;
	//std::cout << "Msub = " << remove_column(M,0) << std::endl;
	//std::cout << ublas::subrange(M,0,2,2,4) << std::endl;

	mat Q, R;
	bool bOk = qr(M, Q, R);
	std::cout << "OK = " << bOk << std::endl;

	std::cout << "\nQ = " << Q << std::endl;
	std::cout << "det(Q) = " << determinant(Q) << std::endl;
	std::cout << "R = " << R << std::endl;

	std::cout << "\nQR = " << ublas::prod(Q,R) << std::endl;
	std::cout << "Q'Q = " << ublas::prod(ublas::trans(Q),Q) << std::endl;
	return 0;
}
