/**
 * advanced linalg helpers (which depend on lapack/e)
 * @author: Tobias Weber <tobias.weber@tum.de>
 * @date: 30-apr-2013 - 2018
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

#include "linalg2.h"
#include "linalg2_impl.h"

#ifndef NO_LAPACK

namespace tl {

template bool qr(const ublas::matrix<double>& M,
	ublas::matrix<double>& Q, ublas::matrix<double>& R);


template bool eigenvec(const ublas::matrix<double>& mat,
	std::vector<ublas::vector<double> >& evecs_real, std::vector<ublas::vector<double>>& evecs_imag,
	std::vector<double>& evals_real, std::vector<double>& evals_imag);

template bool eigenvec_cplx(const ublas::matrix<std::complex<double>>& mat,
	std::vector<ublas::vector<std::complex<double>> >& evecs,
	std::vector<std::complex<double>>& evals);


template bool eigenvec_sym(const ublas::matrix<double>& mat,
	std::vector<ublas::vector<double>>& evecs, std::vector<double>& evals);

template bool eigenvec_approxsym(const ublas::matrix<double>& mat,
	std::vector<ublas::vector<double>>& evecs, std::vector<double>& evals);

template bool eigenvec_herm(const ublas::matrix<std::complex<double>>& mat,
	std::vector<ublas::vector<std::complex<double>>>& evecs,
	std::vector<double>& evals);

}

#endif
