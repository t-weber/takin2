/**
 * Wrapper for clipper spacegroups
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date oct-2015
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

#include "spacegroup_clp.h"
#include <ctype.h>
#include "tlibs/log/log.h"

namespace xtl {


SpaceGroupClp::SpaceGroupClp(unsigned int iNum)
	: m_psg(new clipper::Spacegroup(clipper::Spgr_descr(iNum)))
{
	m_strName = m_psg->symbol_hm();
	convert_hm_symbol(m_strName);

	m_strLaue = m_psg->symbol_laue();
	m_crystalsys = get_crystal_system_from_laue_group(m_strLaue.c_str());
}

SpaceGroupClp::~SpaceGroupClp()
{
	if(m_psg) { delete m_psg; m_psg = nullptr; }
}

SpaceGroupClp::SpaceGroupClp(const SpaceGroupClp& sg)
	: m_psg(new clipper::Spacegroup(sg.m_psg->descr())),
	m_strName(sg.GetName()), m_strLaue(sg.GetLaueGroup()),
	m_crystalsys(sg.GetCrystalSystem())
{}

SpaceGroupClp::SpaceGroupClp(SpaceGroupClp&& sg)
{
	m_psg = sg.m_psg;
	sg.m_psg = nullptr;

	m_strName = std::move(sg.m_strName);
	m_strLaue = std::move(sg.m_strLaue);
	m_crystalsys = std::move(sg.m_crystalsys);
}

// direct calculation of allowed bragg peaks
bool SpaceGroupClp::HasReflection(int h, int k, int l) const
{
	if(!m_psg) return false;

	std::vector<ublas::matrix<double>> vecTrafos;
	GetSymTrafos(vecTrafos);
	return is_reflection_allowed<std::vector,
		ublas::matrix<double>, ublas::vector<double>>(h,k,l, vecTrafos).first;
}

// using clipper
bool SpaceGroupClp::HasReflection2(int h, int k, int l) const
{
	if(!m_psg) return false;

	clipper::HKL_class hkl = m_psg->hkl_class(clipper::HKL(h,k,l));
	return !hkl.sys_abs();
}

void SpaceGroupClp::GetSymTrafos(std::vector<ublas::matrix<double>>& vecTrafos) const
{
	get_symtrafos<clipper::ftype>(*m_psg, vecTrafos, SymTrafoType::ALL);
}
void SpaceGroupClp::GetInvertingSymTrafos(std::vector<ublas::matrix<double>>& vecTrafos) const
{
	get_symtrafos<clipper::ftype>(*m_psg, vecTrafos, SymTrafoType::INVERTING);
}
void SpaceGroupClp::GetPrimitiveSymTrafos(std::vector<ublas::matrix<double>>& vecTrafos) const
{
	get_symtrafos<clipper::ftype>(*m_psg, vecTrafos, SymTrafoType::PRIMITIVE);
}
void SpaceGroupClp::GetCenteringSymTrafos(std::vector<ublas::matrix<double>>& vecTrafos) const
{
	get_symtrafos<clipper::ftype>(*m_psg, vecTrafos, SymTrafoType::CENTERING);
}


// -----------------------------------------------------------------------------


static t_mapSpaceGroups g_mapSpaceGroups;

void init_space_groups()
{
	if(!g_mapSpaceGroups.empty())
	{
		tl::log_warn("Space Groups have already been initialised.");
		return;
	}

	typedef t_mapSpaceGroups::value_type t_val;

	for(int iSg=1; iSg<=230; ++iSg)
	{
		SpaceGroupClp sg(iSg);
		g_mapSpaceGroups.insert(t_val(sg.GetName(), std::move(sg)));
	}
}

const t_mapSpaceGroups* get_space_groups()
{
	if(g_mapSpaceGroups.empty())
	{
		tl::log_warn("Space Groups had not been initialised properly.");
		init_space_groups();
	}

	return &g_mapSpaceGroups;
}

}
