/**
 * boost-specific hacks
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jun-2017
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_BOOST_HACKS_H__
#define __TLIBS_BOOST_HACKS_H__


/**
 * include declarations that were moved out of array.hpp
 */
#include <boost/version.hpp>
#if BOOST_VERSION >= 106400
	#include <boost/serialization/array_wrapper.hpp>
#endif


#endif
