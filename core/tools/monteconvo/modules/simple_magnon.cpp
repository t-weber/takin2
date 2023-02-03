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

#include "simple_magnon.h"
#include "tlibs/string/string.h"
#include "tlibs/log/log.h"
#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"
#include "tlibs/phys/neutrons.h"
#include <fstream>
#include <list>

using t_real = t_real_reso;


t_real SqwMagnon::ferro_disp(t_real dq, t_real dD, t_real doffs)
{
	return dq*dq * dD + doffs;
}


t_real SqwMagnon::antiferro_disp(t_real dq, t_real dD, t_real doffs)
{
	return std::abs(dq)*dD + doffs;
}


SqwMagnon::SqwMagnon(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr)
		tl::log_warn("Cannot open magnon config file \"", pcFile, "\".");

	// if a file is given, load the parameters
	if(!!ifstr)
	{
		std::string strLine;
		while(std::getline(ifstr, strLine))
		{
			std::vector<std::string> vecToks;
			tl::get_tokens<std::string>(strLine, std::string("=,"), vecToks);
			std::for_each(vecToks.begin(), vecToks.end(), [](std::string& str) {tl::trim(str); });

			if(vecToks.size() == 0) continue;

			if(vecToks[0] == "G") m_vecBragg = tl::make_vec({tl::str_to_var_parse<t_real>(vecToks[1]), tl::str_to_var_parse<t_real>(vecToks[2]), tl::str_to_var_parse<t_real>(vecToks[3])});
			else if(vecToks[0] == "disp") m_iWhichDisp = tl::str_to_var<decltype(m_iWhichDisp)>(vecToks[1]);

			else if(vecToks[0] == "D") m_dD = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "offs") m_dOffs = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "E_HWHM") m_dE_HWHM = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "S0") m_dS0 = tl::str_to_var_parse<t_real>(vecToks[1]);

			else if(vecToks[0] == "inc_amp") m_dIncAmp = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "inc_sig") m_dIncSig = tl::str_to_var_parse<t_real>(vecToks[1]);

			else if(vecToks[0] == "T") m_dT = tl::str_to_var_parse<t_real>(vecToks[1]);
		}
	}

	// just use 100 if nothing else is defined
	if(m_vecBragg.size() < 3)
		m_vecBragg = tl::make_vec({1., 0., 0.});

	m_bOk = 1;
}


/**
 * dispersion E(Q)
 */
std::tuple<std::vector<t_real>, std::vector<t_real>>
SqwMagnon::disp(t_real dh, t_real dk, t_real dl) const
{
	dh -= m_vecBragg[0];
	dk -= m_vecBragg[1];
	dl -= m_vecBragg[2];

	t_real (*pDisp)(t_real, t_real, t_real) = nullptr;
	switch(m_iWhichDisp)
	{
		case 0: pDisp = &ferro_disp; break;
		case 1: pDisp = &antiferro_disp; break;
	}

	if(!pDisp)
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());

	t_real dq = std::sqrt(dh*dh + dk*dk + dl*dl);
	t_real dE = pDisp(dq, m_dD, m_dOffs);
	t_real dW = t_real(1);

	return std::make_tuple(std::vector<t_real>({dE, -dE}),
		std::vector<t_real>({dW, dW}));
}


/**
 * dynamical structure factor S(Q,E)
 */
t_real SqwMagnon::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	std::vector<t_real> vecE0, vecW;
	std::tie(vecE0, vecW) = disp(dh, dk, dl);

	t_real dInc = 0.;
	if(!tl::float_equal<t_real>(m_dIncAmp, 0.))
		dInc = tl::gauss_model<t_real>(dE, 0., m_dIncSig, m_dIncAmp, 0.);

	t_real dS = 0;
	if(vecE0.size())
	{
		for(std::size_t i=0; i<vecE0.size(); ++i)
			dS += std::abs(tl::DHO_model<t_real>(dE, m_dT, vecE0[i], m_dE_HWHM, vecW[i], 0.));
		dS *= m_dS0;
	}

	return dS + dInc;
}


std::vector<SqwBase::t_var> SqwMagnon::GetVars() const
{
	std::vector<SqwBase::t_var> vecVars;

	vecVars.push_back(SqwBase::t_var{"G", "vector", vec_to_str(m_vecBragg)});
	vecVars.push_back(SqwBase::t_var{"disp", "uint", tl::var_to_str(m_iWhichDisp)});

	vecVars.push_back(SqwBase::t_var{"D", "real", tl::var_to_str(m_dD)});
	vecVars.push_back(SqwBase::t_var{"offs", "real", tl::var_to_str(m_dOffs)});
	vecVars.push_back(SqwBase::t_var{"E_HWHM", "real", tl::var_to_str(m_dE_HWHM)});
	vecVars.push_back(SqwBase::t_var{"S0", "real", tl::var_to_str(m_dS0)});

	vecVars.push_back(SqwBase::t_var{"inc_amp", "real", tl::var_to_str(m_dIncAmp)});
	vecVars.push_back(SqwBase::t_var{"inc_sig", "real", tl::var_to_str(m_dIncSig)});

	vecVars.push_back(SqwBase::t_var{"T", "real", tl::var_to_str(m_dT)});

	return vecVars;
}


void SqwMagnon::SetVars(const std::vector<SqwBase::t_var>& vecVars)
{
	if(vecVars.size() == 0)
		return;

	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "G") m_vecBragg = str_to_vec<decltype(m_vecBragg)>(strVal);
		else if(strVar == "disp") m_iWhichDisp = tl::str_to_var<decltype(m_iWhichDisp)>(strVal);

		else if(strVar == "D") m_dD = tl::str_to_var<decltype(m_dD)>(strVal);
		else if(strVar == "offs") m_dOffs = tl::str_to_var<decltype(m_dOffs)>(strVal);
		else if(strVar == "E_HWHM") m_dE_HWHM = tl::str_to_var<decltype(m_dE_HWHM)>(strVal);
		else if(strVar == "S0") m_dS0 = tl::str_to_var<decltype(m_dS0)>(strVal);

		else if(strVar == "inc_amp") m_dIncAmp = tl::str_to_var<decltype(m_dIncAmp)>(strVal);
		else if(strVar == "inc_sig") m_dIncSig = tl::str_to_var<decltype(m_dIncSig)>(strVal);

		else if(strVar == "T") m_dT = tl::str_to_var<decltype(m_dT)>(strVal);
	}
}


SqwBase* SqwMagnon::shallow_copy() const
{
	SqwMagnon *pCpy = new SqwMagnon();
	*static_cast<SqwBase*>(pCpy) = *static_cast<const SqwBase*>(this);

	pCpy->m_vecBragg = m_vecBragg;
	pCpy->m_iWhichDisp = m_iWhichDisp;

	pCpy->m_dD = m_dD;
	pCpy->m_dOffs = m_dOffs;
	pCpy->m_dS0 = m_dS0;
	pCpy->m_dE_HWHM = m_dE_HWHM;

	pCpy->m_dIncAmp = m_dIncAmp;
	pCpy->m_dIncSig = m_dIncSig;

	pCpy->m_dT = m_dT;

	return pCpy;
}
