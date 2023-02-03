/**
 * libcrystal main interface
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2017
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

#define _FF_NO_SINGLETON
#define _SGR_NO_SINGLETON
#include "libs/spacegroups/spacegroup.h"
#include "libs/formfactors/formfact.h"

#include "libxtl.h"

using t_real = double;

xtl::SpaceGroups<t_real>* g_pSGs = nullptr;
xtl::ScatlenList<t_real>* g_pSLs = nullptr;


/**
 * loads needed tables
 */
int xtl_init()
{
	if(g_pSGs || g_pSLs)
		xtl_deinit();

	g_pSGs = new xtl::SpaceGroups<t_real>("libxtl.xml", "xtl");
	g_pSLs = new xtl::ScatlenList<t_real>("libxtl.xml", "xtl");

	return 0;
}


/**
 * unloads tables
 */
void xtl_deinit()
{
	if(g_pSGs)
	{
		delete g_pSGs;
		g_pSGs = nullptr;
	}

	if(g_pSLs)
	{
		delete g_pSLs;
		g_pSLs = nullptr;
	}
}
