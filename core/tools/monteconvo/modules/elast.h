/**
 * pre-defined sqw modules
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015 -- 2018
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

#ifndef __MCONV_SQWMOD_ELAST_H__
#define __MCONV_SQWMOD_ELAST_H__

#include <list>
#include <unordered_map>

#include "tlibs/helper/boost_hacks.h"
#include <boost/numeric/ublas/vector.hpp>

#include "tlibs/math/math.h"
#include "tlibs/math/kd.h"
#include "tlibs/file/loaddat.h"
#include "../../res/defs.h"
#include "../sqwbase.h"


namespace ublas = boost::numeric::ublas;


struct ElastPeak
{
	t_real_reso h, k, l;
	t_real_reso dSigQ, dSigE;
	t_real_reso dS;
};


/**
 * Bragg peaks
 */
class SqwElast : public SqwBase
{
protected:
	bool m_bLoadedFromFile = false;
	std::list<ElastPeak> m_lstPeaks;

public:
	SqwElast() { SqwBase::m_bOk = true; }
	SqwElast(const char* pcFile);
	virtual t_real_reso operator()(t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override;

	void AddPeak(t_real_reso h, t_real_reso k, t_real_reso l, t_real_reso dSigQ, t_real_reso dSigE, t_real_reso dS);

	virtual std::vector<SqwBase::t_var> GetVars() const override;
	virtual void SetVars(const std::vector<SqwBase::t_var>&) override;

	virtual SqwBase* shallow_copy() const override;
};


#endif
