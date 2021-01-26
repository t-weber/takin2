/**
 * libcrystal main interface
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2017
 * @license GPLv2
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
