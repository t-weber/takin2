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

#include "kdtree.h"
#include "tlibs/string/string.h"
#include "tlibs/log/log.h"
#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"
#include "tlibs/phys/neutrons.h"
#include <fstream>
#include <list>

using t_real = t_real_reso;


SqwKdTree::SqwKdTree(const char* pcFile)
{
	if(pcFile)
		m_bOk = open(pcFile);
}


bool SqwKdTree::open(const char* pcFile)
{
	m_kd = std::make_shared<tl::Kd<t_real>>();

	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

	std::list<std::vector<t_real>> lstPoints;
	std::size_t iCurPoint = 0;
	while(!ifstr.eof())
	{
		std::string strLine;
		std::getline(ifstr, strLine);
		tl::trim(strLine);

		if(strLine.length() == 0)
			continue;

		if(strLine[0] == '#')
		{
			strLine[0] = ' ';
			m_mapParams.insert(tl::split_first(strLine, std::string(":"), 1));
			continue;
		}

		std::vector<t_real> vecSqw;
		tl::get_tokens<t_real>(strLine, std::string(" \t"), vecSqw);
		if(vecSqw.size() != 5)
		{
			tl::log_err("Need h,k,l,E,S data.");
			return false;
		}

		lstPoints.push_back(vecSqw);
		++iCurPoint;
	}

	tl::log_info("Loaded ",  iCurPoint, " S(q,w) points.");
	m_kd->Load(lstPoints, 4);
	tl::log_info("Generated k-d tree.");

	return true;
}


t_real SqwKdTree::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	// meV and rlu units will have equal scaling in the kd tree!
	std::vector<t_real> vechklE = {dh, dk, dl, dE};

	// Warning: Will return 0 when bounding box has one of the four dimensions is 0
	if(!m_kd->IsPointInGrid(vechklE))
	{
		//std::cout << "not in grid: (" << dh << " " << dk << " " << dl << " " << dE << ")" << std::endl;
		return 0.;
	}

	std::vector<t_real> vec = m_kd->GetNearestNode(vechklE);
	//std::cout << "querying (" << dh << " " << dk << " " << dl << " " << dE << "), result: "
	//	<< vec[0] << " " << vec[1] << " " << vec[2] << " " << vec[3] << std::endl;

	return vec[4];
}


std::vector<SqwBase::t_var> SqwKdTree::GetVars() const
{
	std::vector<SqwBase::t_var> vecVars;

	return vecVars;
}


void SqwKdTree::SetVars(const std::vector<SqwBase::t_var>&)
{
}


SqwBase* SqwKdTree::shallow_copy() const
{
	SqwKdTree *pTree = new SqwKdTree();
	*static_cast<SqwBase*>(pTree) = *static_cast<const SqwBase*>(this);

	pTree->m_mapParams = m_mapParams;
	pTree->m_kd = m_kd;

	return pTree;
}
