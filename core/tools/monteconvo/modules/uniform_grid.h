/**
 * S(q,w) module for a pre-calculated grid (file format version 2)
 * @author Tobias Weber <tweber@ill.fr>
 * @date 06-jan-2020
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

#ifndef __MCONV_SQW_GRID_VER2_H__
#define __MCONV_SQW_GRID_VER2_H__

#include "tools/monteconvo/sqwbase.h"
#include "tlibs/math/linalg.h"


class SqwUniformGrid : public SqwBase
{
	public:
		using SqwBase::t_var;
		using t_real = t_real_reso;
		using t_vec = tl::ublas::vector<t_real>;

	protected:
		// temperature for Bose factor
		t_real m_dT = t_real(100);

		// Bose cutoff
		t_real m_dcut = t_real(0.02);

		// peak width for creation and annihilation
		t_real m_dSigma = t_real(0.05);

		// S(q,E) scaling factor
		t_real m_dS0 = t_real(1.);

		// incoherent amplitude and width
		t_real m_dIncAmp = t_real(0.);
		t_real m_dIncSigma = t_real(0.05);

		// grid data file
		std::string m_strDataFile;
		std::size_t m_indexBlockOffset = 0;

		t_real m_hmin=0., m_hmax=0., m_hstep=0.;
		t_real m_kmin=0., m_kmax=0., m_kstep=0.;
		t_real m_lmin=0., m_lmax=0., m_lstep=0.;
		std::size_t m_hsize=0, m_ksize=0, m_lsize=0;


	public:
		SqwUniformGrid();
		SqwUniformGrid(const std::string& strDatFile);
		virtual ~SqwUniformGrid();

		virtual std::tuple<std::vector<t_real>, std::vector<t_real>>
			disp(t_real dh, t_real dk, t_real dl) const override;
		virtual t_real operator()(t_real dh, t_real dk, t_real dl, t_real dE) const override;

		virtual std::vector<t_var> GetVars() const override;
		virtual void SetVars(const std::vector<t_var>&) override;
		virtual bool SetVarIfAvail(const std::string& strKey, const std::string& strNewVal) override;

		virtual SqwBase* shallow_copy() const override;
};

#endif
