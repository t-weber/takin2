/**
 * globals
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 20-mar-2015
 * @license GPLv2
 */

#ifndef __XTL_GLOBALS_H__
#define __XTL_GLOBALS_H__

#include <vector>
#include <string>


//using t_real_glob = float;
using t_real_glob = double;


extern unsigned int g_iPrec;
extern unsigned int g_iPrecGfx;

extern t_real_glob g_dEps;
extern t_real_glob g_dEpsGfx;

extern unsigned int g_iMaxThreads;

extern void add_resource_path(const std::string& strPath, bool bToBack=1);
extern std::string find_resource(const std::string& strFile, bool bLogErr=1);
extern std::vector<std::string> find_resource_dirs(const std::string& strDir, bool bLogErr=1);

extern unsigned int get_max_threads();


#endif
