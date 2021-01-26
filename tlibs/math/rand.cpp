/**
 * random numbers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 16-aug-2013
 * @license GPLv2 or GPLv3
 */

#include "rand.h"
#include "../log/log.h"

#include <cstdlib>
#include <exception>
#include <time.h>
#include <sys/time.h>

namespace tl {
// thread-global engine to use if a thread-local has not been seeded
static std::mt19937 g_randeng_fallback;
static bool g_bFallbackInited = 0;

static thread_local std::mt19937/*_64*/ g_randeng;
static thread_local bool g_bHasEntropy = 0;
static thread_local bool g_bIsSeeded = 0;


std::mt19937& get_randeng()
{
	if(g_bIsSeeded)
		return g_randeng;
	return g_randeng_fallback;
}


unsigned int get_rand_seed()
{
	// seed 0: random device
	unsigned int uiSeed0 = 0;
	try
	{
		std::random_device rnd;
		g_bHasEntropy = (rnd.entropy() != 0);
		uiSeed0 = rnd();
	}
	catch(const std::exception& ex)
	{
		log_debug(ex.what());
		uiSeed0 = 0;
	}

	// seed 1: time based
	struct timeval timev;
	gettimeofday(&timev, 0);
	unsigned int uiSeed1 = timev.tv_sec ^ timev.tv_usec;

	// total seed
	unsigned int uiSeed = uiSeed0 ^ uiSeed1;
	return uiSeed;
}


void init_rand()
{
	init_rand_seed(get_rand_seed());
}


void init_rand_seed(unsigned int uiSeed)
{
	std::string strEntr;
	if(!g_bHasEntropy)
		strEntr = ", but entropy is zero";
	log_debug("Random seed: ", uiSeed, strEntr, ".");

	srand(uiSeed);
	g_randeng = std::mt19937/*_64*/(uiSeed);
	g_bIsSeeded = 1;

	// copy first engine to thread-global fallback
	if(!g_bFallbackInited)
	{
		g_randeng_fallback = g_randeng;
		g_bFallbackInited = 1;
	}
}


unsigned int simple_rand(unsigned int iMax)
{
	return rand() % iMax;
}

}
