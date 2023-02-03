// gcc -o domains domains.cpp -std=c++11 -lstdc++ -lm -I../.. -I/usr/include/QtGui/ ../../libs/spacegroups/spacegroup_clp.cpp ../../libs/spacegroups/crystalsys.cpp -DNO_QT -lclipper-core ../../tlibs/log/log.cpp

/**
 * generates positions based on space group
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
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

#include <clipper/clipper.h>
#include <vector>
#include <sstream>
#include "tlibs/math/linalg.h"
#include "tlibs/math/linalg_ops.h"
#include "tlibs/string/string.h"
#include "libs/spacegroups/spacegroup_clp.h"


typedef tl::ublas::vector<double> t_vec;
typedef tl::ublas::matrix<double> t_mat;

void gen_dirs()
{
	std::string strSg;
	std::cout << "Enter spacegroup: ";
	std::getline(std::cin, strSg);
	tl::trim(strSg);
	clipper::Spgr_descr dsc(strSg);

	const int iSGNum = dsc.spacegroup_number();
	if(iSGNum <= 0)
	{
		std::cerr << "Error: Unknown spacegroup." << std::endl;
		return;
	}
	std::cout << "Spacegroup number: " << iSGNum << std::endl;


	clipper::Spacegroup sg(dsc);
	std::vector<t_mat> vecTrafos;
	get_symtrafos(sg, vecTrafos);
	std::cout << vecTrafos.size() << " symmetry operations." << std::endl;


	t_vec vecDir(4);
	std::cout << "Enter direction: ";
	std::cin >> vecDir[0] >> vecDir[1] >> vecDir[2];
	vecDir[3] = 0.;		// no translations, only point groups


	std::vector<t_vec> vecNewDirs;
	std::cout << "\nall transformations:" << std::endl;
	for(const t_mat& matTrafo : vecTrafos)
	{
		t_vec vecDirNew = ublas::prod(matTrafo, vecDir);
		std::cout << vecDirNew << " (from trafo " << matTrafo << ")" << std::endl;

		bool bHasDir = 0;
		for(const t_vec& vec : vecNewDirs)
			if(tl::vec_equal(vec, vecDirNew))
			{
				bHasDir = 1;
				break;
			}
		if(!bHasDir)
			vecNewDirs.push_back(vecDirNew);
	}

	std::cout << "\nunique transformations:" << std::endl;
	for(const t_vec& vec : vecNewDirs)
	{
		std::cout << vec << std::endl;
	}
}


int main()
{
	try
	{
		gen_dirs();
	}
	catch(const clipper::Message_fatal& err)
	{
		std::cerr << "Error in spacegroup." << std::endl;
	}
	return 0;
}
