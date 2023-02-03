/**
 * monte carlo convolution tool
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

#include "elast.h"
#include "tlibs/string/string.h"
#include "tlibs/log/log.h"
#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"
#include "tlibs/phys/neutrons.h"
#include <fstream>
#include <list>

using t_real = t_real_reso;


SqwElast::SqwElast(const char* pcFile) : m_bLoadedFromFile(true)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr)
	{
		m_bLoadedFromFile = false;
		tl::log_err("Cannot open config file.");
		return;
	}

	std::string strLine;
	while(std::getline(ifstr, strLine))
	{
		tl::trim(strLine);
		if(strLine.length()==0 || strLine[0]=='#')
			continue;

		std::istringstream istr(strLine);
		t_real h=0., k=0. ,l=0., dSigQ=0., dSigE=0., dS=0.;
		istr >> h >> k >> l >> dSigQ >> dSigE >> dS;

		AddPeak(h,k,l, dSigQ, dSigE, dS);
	}

	tl::log_info("Number of elastic peaks: ", m_lstPeaks.size());
	SqwBase::m_bOk = true;
}

void SqwElast::AddPeak(t_real h, t_real k, t_real l, t_real dSigQ, t_real dSigE, t_real dS)
{
	ElastPeak pk;
	pk.h = h; pk.k = k; pk.l = l;
	pk.dSigQ = dSigQ; pk.dSigE = dSigE;
	pk.dS = dS;
	m_lstPeaks.push_back(std::move(pk));
}

t_real SqwElast::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	const ublas::vector<t_real> vecCur = tl::make_vec({dh, dk, dl});

	if(!m_bLoadedFromFile)	// use nearest integer bragg peak
	{
		const ublas::vector<t_real> vecPt = tl::make_vec({std::round(dh), std::round(dk), std::round(dl)});

		const t_real dDistQ = ublas::norm_2(vecPt-vecCur);
		const t_real dSigmaQ = 0.02;
		const t_real dSigmaE = 0.02;

		return tl::gauss_model<t_real>(dDistQ, 0., dSigmaQ, 1., 0.) *
			tl::gauss_model<t_real>(dE, 0., dSigmaE, 1., 0.);
	}
	else	// use bragg peaks from config file
	{
		t_real dS = 0.;

		for(const ElastPeak& pk : m_lstPeaks)
		{
			const ublas::vector<t_real> vecPk = tl::make_vec({pk.h, pk.k, pk.l});
			const t_real dDistQ = ublas::norm_2(vecPk-vecCur);

			dS += pk.dS * tl::gauss_model<t_real>(dDistQ, 0., pk.dSigQ, 1., 0.) *
				tl::gauss_model<t_real>(dE, 0., pk.dSigE, 1., 0.);
		}

		return dS;
	}
}


std::vector<SqwBase::t_var> SqwElast::GetVars() const
{
	std::vector<SqwBase::t_var> vecVars;
	return vecVars;
}


void SqwElast::SetVars(const std::vector<SqwBase::t_var>&)
{
}


SqwBase* SqwElast::shallow_copy() const
{
	SqwElast *pElast = new SqwElast();
	*static_cast<SqwBase*>(pElast) = *static_cast<const SqwBase*>(this);

	pElast->m_bLoadedFromFile = m_bLoadedFromFile;
	pElast->m_lstPeaks = m_lstPeaks;	// not a shallow copy!
	return pElast;
}
