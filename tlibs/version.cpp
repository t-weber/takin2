/**
 * tlibs
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2012-2020
 * @license GPLv2 or GPLv3
 */

#include "version.h"
#include <string>

namespace tl {


const char* get_tlibs_version()
{
	return TLIBS_VERSION;
}

const char* get_tlibs_infos()
{
	return "This is the tlibs template library.\n"
		"Written by Tobias Weber <tobias.weber@tum.de>, 2012 - 2021.\n"
		"License: GPLv2 or GPLv3.";
}

/**
 * Check if supplied string matches with compiled one
 * (to check if .so library and header files match)
 */
bool check_tlibs_version(const char* pcHdrVer)
{
	return std::string(TLIBS_VERSION) == std::string(pcHdrVer);
}


}
