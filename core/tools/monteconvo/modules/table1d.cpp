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

#include "table1d.h"
#include "tlibs/string/string.h"
#include "tlibs/log/log.h"
#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"
#include "tlibs/phys/neutrons.h"
#include <fstream>
#include <list>

using t_real = t_real_reso;


SqwTable1d::SqwTable1d(const char* pcFile)
{
	m_bOk = false;
	if(pcFile)
		m_bOk = open(pcFile);
}


bool SqwTable1d::open(const char* pcFile)
{
	tl::log_debug("Loading \"", pcFile, "\"", ".");
	m_dat = std::make_shared<tl::DatFile<t_real>>(pcFile);
	m_bOk = m_dat->IsOk();
	CreateKd();
	return m_bOk;
}


void SqwTable1d::CreateKd()
{
	if(!m_bOk)
	{
		tl::log_err("Data table not yet loaded.");
		return;
	}

	m_kd = std::make_shared<tl::Kd<t_real>>();
	std::list<std::vector<t_real>> lstPoints;

	t_real minq = std::numeric_limits<t_real>::max();
	t_real maxq = -std::numeric_limits<t_real>::max();
	t_real minE = std::numeric_limits<t_real>::max();
	t_real maxE = -std::numeric_limits<t_real>::max();

	for(std::size_t iRow=0; iRow<m_dat->GetRowCount(); ++iRow)
	{
		t_real q = m_dat->GetColumn(m_qcol)[iRow];
		t_real E = m_dat->GetColumn(m_Ecol)[iRow];
		t_real S = m_dat->GetColumn(m_Scol)[iRow];

		minq = std::min(minq, q);
		maxq = std::max(maxq, q);
		minE = std::min(minE, E);
		maxE = std::max(maxE, E);

		lstPoints.emplace_back(std::vector<t_real>{{ q, E, S }});
	}

	tl::log_info("Loaded ", m_dat->GetRowCount(), " S(q,w) points.");
	tl::log_info("q range: ", minq, "..", maxq, ", E range: ", minE, "..", maxE, ".");

	m_kd->Load(lstPoints, 2);
	tl::log_info("Generated k-d tree.");
}


t_real SqwTable1d::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	if(!m_bOk)
		return 0.;

	// get reduced q
	dh -= m_G[0]; dk -= m_G[1]; dl -= m_G[2];
	t_real dq = std::sqrt(dh*dh + dk*dk + dl*dl);

	// meV and rlu units will have equal scaling in the kd tree!
	std::vector<t_real> vecqE{{dq, dE}};

	if(!m_kd->IsPointInGrid(vecqE))
		return 0.;

	std::vector<t_real> vec = m_kd->GetNearestNode(vecqE);
	return vec[2];
}


std::vector<SqwBase::t_var> SqwTable1d::GetVars() const
{
	std::ostringstream ostr;
	ostr.precision(g_iPrec);
	ostr << m_G[0] << " " << m_G[1] << " " << m_G[2];


	std::vector<SqwBase::t_var> vecVars;

	vecVars.push_back(SqwBase::t_var{"q_column", "uint", tl::var_to_str(m_qcol)});
	vecVars.push_back(SqwBase::t_var{"E_column", "uint", tl::var_to_str(m_Ecol)});
	vecVars.push_back(SqwBase::t_var{"S_column", "uint", tl::var_to_str(m_Scol)});
	vecVars.push_back(SqwBase::t_var{"G", "vector", ostr.str()});

	return vecVars;
}


void SqwTable1d::SetVars(const std::vector<SqwBase::t_var>& vecVars)
{
	if(vecVars.size() == 0)
		return;

	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "q_column") m_qcol = tl::str_to_var<decltype(m_qcol)>(strVal);
		else if(strVar == "E_column") m_Ecol = tl::str_to_var<decltype(m_Ecol)>(strVal);
		else if(strVar == "S_column") m_Scol = tl::str_to_var<decltype(m_Scol)>(strVal);
		else if(strVar == "G")
		{
			std::istringstream istr(strVal);
			istr >> m_G[0] >> m_G[1] >> m_G[2];
		}
	}

	CreateKd();
}


SqwBase* SqwTable1d::shallow_copy() const
{
	SqwTable1d *pTab = new SqwTable1d();
	*static_cast<SqwBase*>(pTab) = *static_cast<const SqwBase*>(this);

	pTab->m_dat = m_dat;
	pTab->m_kd = m_kd;

	return pTab;
}
