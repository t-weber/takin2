/**
 * tlibs
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2012-2021
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_VER_H__
#define __TLIBS_VER_H__

#define TLIBS_VERSION "0.9.7"

namespace tl {

extern const char* get_tlibs_version();
extern const char* get_tlibs_infos();
extern bool check_tlibs_version(const char* pcHdrVer);

}
#endif
