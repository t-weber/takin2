/**
 * @author Tobias Weber <tobias.weber@tum.de>
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

// gcc -DNO_QT -I. -I../.. -o tst_refl tst_refl.cpp ../../libs/spacegroups/spacegroup.cpp ../../libs/spacegroups/crystalsys.cpp ../../tlibs/log/log.cpp ../../libs/globals.cpp -lstdc++ -std=c++11 -lm -lboost_iostreams -lboost_filesystem -lboost_system

#include <iostream>
#include "libs/spacegroups/spacegroup.h"

using t_real = double;
using t_mapSpaceGroups = xtl::SpaceGroups<t_real>::t_mapSpaceGroups;


void check_allowed_refls()
{
	std::shared_ptr<const xtl::SpaceGroups<t_real>> sgs = xtl::SpaceGroups<t_real>::GetInstance();

	const int HKL_MAX = 10;
	unsigned int iSG = 0;

	const t_mapSpaceGroups *pSGs = sgs->get_space_groups();
	for(const t_mapSpaceGroups::value_type& sg : *pSGs)
	{
		++iSG;
		//if(iSG < 200) continue;

		const std::string& strName = sg.second.GetName();
		std::cout << "Checking (" << iSG << ") " << strName << " ... ";

		for(int i=-HKL_MAX; i<=HKL_MAX; ++i)
		for(int j=-HKL_MAX; j<=HKL_MAX; ++j)
		for(int k=-HKL_MAX; k<=HKL_MAX; ++k)
		{
			//std::cout << "(" << i << j << k << "): " << 
			//	std::boolalpha << sg.second.HasReflection(i,j,k) << std::endl;

			if(sg.second.HasReflection(i,j,k) != sg.second.HasReflection2(i,j,k))
			{
				std::cout << "Failed at " << i << " " << j << " " << k << std::endl;
				return;
			}
		}

		std::cout << "OK" << std::endl;
	}
}

int main()
{
	check_allowed_refls();
	return 0;
}
