/**
 * Compiler- and system-specific stuff
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_COMPILER_FLAGS_H__
#define __TLIBS_COMPILER_FLAGS_H__


namespace tl {

// normal popen is not thread-safe on all systems
void *my_popen(const char* pcCmd, const char* pcType="w");
int my_pclose(void*);

}

#endif
