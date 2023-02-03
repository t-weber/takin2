/**
 * fftw interface
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date August 2012
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

#ifndef __FOURIER_FFTW__
#define __FOURIER_FFTW__

#include "dft.h"

namespace tl {


class FFTw : public Fourier_base<double>
{
	protected:
		std::size_t m_iSize;

		void *m_pIn, *m_pOut;
		void *m_pPlan, *m_pPlan_inv;


	public:
		FFTw(unsigned int iSize);
		virtual ~FFTw();

		virtual void trafo(const double* pInR, const double *pInI,
			double *pOutR, double *pOutI, bool bInv=0) override;
};

}

#endif
