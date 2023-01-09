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

#ifndef __MCONV_SQWMOD_MAGNON_H__
#define __MCONV_SQWMOD_MAGNON_H__

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


/**
 * simple magnon model
 */
class SqwMagnon : public SqwBase
{
private:
	SqwMagnon() {};

protected:
	static t_real_reso ferro_disp(t_real_reso dq, t_real_reso dD, t_real_reso doffs);
	static t_real_reso antiferro_disp(t_real_reso dq, t_real_reso dD, t_real_reso doffs);

protected:
	unsigned short m_iWhichDisp = 0;		// 0: ferro, 1: antiferro
	ublas::vector<t_real_reso> m_vecBragg;

	t_real_reso m_dD = 1., m_dOffs = 0.;
	t_real_reso m_dE_HWHM = 0.1;
	t_real_reso m_dS0 = 1.;

	t_real_reso m_dIncAmp = 0., m_dIncSig = 0.1;
	t_real_reso m_dT = 300.;

public:
	SqwMagnon(const char* pcFile);
	virtual ~SqwMagnon() = default;

	virtual std::tuple<std::vector<t_real_reso>, std::vector<t_real_reso>>
		disp(t_real_reso dh, t_real_reso dk, t_real_reso dl) const override;
	virtual t_real_reso operator()(t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override;

	const ublas::vector<t_real_reso>& GetBragg() const { return m_vecBragg; }

	virtual std::vector<SqwBase::t_var> GetVars() const override;
	virtual void SetVars(const std::vector<SqwBase::t_var>&) override;

	virtual SqwBase* shallow_copy() const override;
};


#endif
