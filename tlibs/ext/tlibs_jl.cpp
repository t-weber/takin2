/**
 * julia module
 * gcc -O2 -march=native -I. -I.. -I/usr/include/julia -I/usr/include/root -std=c++11 -shared -fPIC -o tlibs_jl.so ../log/log.cpp tlibs_jl.cpp -lstdc++ -lboost_system -lboost_iostreams -L/usr/lib64/root -lMinuit2
 * gcc -O2 -march=native -I. -I.. -I/usr/include/julia -I/usr/local/include -std=c++11 -shared -fPIC -o tlibs_jl.so ../log/log.cpp tlibs_jl.cpp -lstdc++ -lboost_system -lboost_iostreams -L/usr/local/lib64 -lMinuit2 -lgomp
 *
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 23-apr-2017
 * @license GPLv2 or GPLv3
 */

#define TLIBS_INC_HDR_IMPLS
#include "version.h"
#include "ext/jl.h"
#include "log/log.h"
#include "file/loadinstr.h"
#include "fit/minuit.h"


using t_real = double;


bool g_bDebug = 0;


extern "C" void load_tlibs(int bDebug)
{
	g_bDebug = (bDebug != 0);
	tl::log_debug.SetEnabled(g_bDebug);

	tl::log_debug("Loaded tlibs version ", TLIBS_VERSION, ".");
}


// ----------------------------------------------------------------------------

/**
 * loads an instrument data file
 */
extern "C" jl_array_t* load_instr(const char* pcFile)
{
	tl::FileInstrBase<t_real>* pInstr = tl::FileInstrBase<t_real>::LoadInstr(pcFile);
	if(!pInstr)
	{
		jl_array_t *pArrNull = jl_alloc_array_1d(
			jl_apply_array_type(reinterpret_cast<jl_value_t*>(jl_any_type), 1), 0);

		tl::log_err("In ", __func__, ": Cannot load ", pcFile, ".");
		return pArrNull;
	}

	// [ column names, data, keys, values ]
	jl_array_t *pArr = jl_alloc_array_1d(
		jl_apply_array_type(reinterpret_cast<jl_value_t*>(jl_any_type), 1), 4);
	jl_array_t** pArrDat = reinterpret_cast<jl_array_t**>(jl_array_data(pArr));

	// data column names
	pArrDat[0] = tl::make_jl_str_arr(pInstr->GetColNames());

	// data matrix
	pArrDat[1] = tl::make_jl_2darr(pInstr->GetData());

	// scan property map
	std::tie(pArrDat[2], pArrDat[3]) = tl::make_jl_strmap_arr(pInstr->GetAllParams());

	return pArr;
}



// ----------------------------------------------------------------------------
// Fitting


/**
 * internal helper function to call fitter using variable args
 */
template<std::size_t iNumArgs>
static inline bool _invoke_fit(void *_pFkt,
	const std::vector<tl::t_real_min>& vecX, const std::vector<tl::t_real_min>& vecY,
	const std::vector<tl::t_real_min>& vecYerr,
	const std::vector<std::string>& vecParamNames,
	std::vector<tl::t_real_min>& vecVals, std::vector<tl::t_real_min>& vecErrs,
	const std::vector<bool>& vecFixed)
{
	auto *pFkt = reinterpret_cast<tl::t_fkt_vararg<t_real, iNumArgs>>(_pFkt);
	return tl::fit<iNumArgs>(pFkt,
		vecX, vecY, vecYerr, vecParamNames, vecVals, vecErrs, &vecFixed, g_bDebug);
}



/**
 * function fitting
 */
extern "C" int fit(void *_pFkt, std::size_t iNumParams,
	const t_real *pX, const t_real *pY, const t_real *pYerr, std::size_t iArrLen,
	jl_array_t *parrParamNames, jl_array_t *parrFixed,
	t_real *pValues, t_real *pErrors)
{
	std::vector<tl::t_real_min> vecX, vecY, vecYerr;

	vecX.reserve(iArrLen);
	vecY.reserve(iArrLen);
	vecYerr.reserve(iArrLen);

	// copy arrays to vectors
	for(std::size_t i=0; i<iArrLen; ++i)
	{
		vecX.push_back(pX[i]);
		vecY.push_back(pY[i]);
		vecYerr.push_back(pYerr[i]);
	}


	// parameter names
	std::vector<std::string> vecParamNames =
		tl::from_jl_arr<std::vector, std::string>(parrParamNames, 1);

	// fixed parameters
	std::vector<std::string> vecFixedParams =
		tl::from_jl_arr<std::vector, std::string>(parrFixed);

	std::vector<bool> vecFixed;
	for(const std::string& strParam : vecParamNames)
	{
		bool bFixed = std::find(vecFixedParams.begin(), vecFixedParams.end(), strParam)
			!= vecFixedParams.end();
		vecFixed.push_back(bFixed);
	}

	// values & errors
	std::vector<tl::t_real_min> vecVals, vecErrs;
	for(std::size_t iIdx=0; iIdx<iNumParams; ++iIdx)
	{
		vecVals.push_back(pValues[iIdx]);
		vecErrs.push_back(pErrors[iIdx]);
	}


	// ------------------------------------------------------------------------
	// fill up missing parameters and hints
	if(vecParamNames.size() < iNumParams)
	{
		for(std::size_t iArg=vecParamNames.size(); iArg<iNumParams; ++iArg)
		{
			std::ostringstream ostrArg;
			ostrArg << "arg_" << iArg;
			vecParamNames.push_back(ostrArg.str());
		}
	}

	while(vecVals.size() < iNumParams) vecVals.push_back(0);
	while(vecErrs.size() < iNumParams) vecErrs.push_back(0);
	while(vecFixed.size() < iNumParams) vecFixed.push_back(0);
	// ------------------------------------------------------------------------


	#define __CALL_FIT(NUM) _invoke_fit<NUM+1>(_pFkt, vecX, vecY, vecYerr, vecParamNames, vecVals, vecErrs, vecFixed)
	bool bOk = 0;

	// stating all needed specialisations of the fit template function
	if(iNumParams == 1) bOk = __CALL_FIT(1);
	else if(iNumParams == 2) bOk = __CALL_FIT(2);
	else if(iNumParams == 3) bOk = __CALL_FIT(3);
	else if(iNumParams == 4) bOk = __CALL_FIT(4);
	else if(iNumParams == 5) bOk = __CALL_FIT(5);
	else if(iNumParams == 6) bOk = __CALL_FIT(6);
	else if(iNumParams == 7) bOk = __CALL_FIT(7);
	else if(iNumParams == 8) bOk = __CALL_FIT(8);
	else if(iNumParams == 9) bOk = __CALL_FIT(9);
	else if(iNumParams == 10) bOk = __CALL_FIT(10);
	else if(iNumParams == 11) bOk = __CALL_FIT(11);
	else if(iNumParams == 12) bOk = __CALL_FIT(12);
	else if(iNumParams == 13) bOk = __CALL_FIT(13);
	else if(iNumParams == 14) bOk = __CALL_FIT(14);
	else if(iNumParams == 15) bOk = __CALL_FIT(15);
	else tl::log_err("In ", __func__, ": Invalid number of function arguments.");


	// write back fitted values & errors
	for(std::size_t iParam=0; iParam<iNumParams; ++iParam)
	{
		pValues[iParam] = vecVals[iParam];
		pErrors[iParam] = vecErrs[iParam];
	}


	return int(bOk);
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Minimisation


/**
 * internal helper function to call fitter using variable args
 */
template<std::size_t iNumArgs>
static inline bool _invoke_minimise(void *_pFkt,
	const std::vector<std::string>& vecParamNames,
	std::vector<tl::t_real_min>& vecVals, std::vector<tl::t_real_min>& vecErrs,
	const std::vector<bool>& vecFixed)
{
	auto *pFkt = reinterpret_cast<tl::t_fkt_vararg<t_real, iNumArgs>>(_pFkt);
	return tl::minimise<iNumArgs>(pFkt, vecParamNames, vecVals, vecErrs, &vecFixed, g_bDebug);
}



/**
 * function minimisation
 */
extern "C" int minimise(void *_pFkt, std::size_t iNumParams,
	jl_array_t *parrParamNames, jl_array_t *parrFixed,
	t_real *pValues, t_real *pErrors)
{
	// parameter names
	std::vector<std::string> vecParamNames =
		tl::from_jl_arr<std::vector, std::string>(parrParamNames);

	// fixed parameters
	std::vector<std::string> vecFixedParams =
		tl::from_jl_arr<std::vector, std::string>(parrFixed);

	std::vector<bool> vecFixed;
	for(const std::string& strParam : vecParamNames)
	{
		bool bFixed = std::find(vecFixedParams.begin(), vecFixedParams.end(), strParam)
			!= vecFixedParams.end();
		vecFixed.push_back(bFixed);
	}

	// values & errors
	std::vector<tl::t_real_min> vecVals, vecErrs;
	for(std::size_t iIdx=0; iIdx<iNumParams; ++iIdx)
	{
		vecVals.push_back(pValues[iIdx]);
		vecErrs.push_back(pErrors[iIdx]);
	}


	// ------------------------------------------------------------------------
	// fill up missing parameters and hints
	if(vecParamNames.size() < iNumParams)
	{
		for(std::size_t iArg=vecParamNames.size(); iArg<iNumParams; ++iArg)
		{
			std::ostringstream ostrArg;
			ostrArg << "arg_" << iArg;
			vecParamNames.push_back(ostrArg.str());
		}
	}

	while(vecVals.size() < iNumParams) vecVals.push_back(0);
	while(vecErrs.size() < iNumParams) vecErrs.push_back(0);
	while(vecFixed.size() < iNumParams) vecFixed.push_back(0);
	// ------------------------------------------------------------------------


	#define __CALL_MINI(NUM) _invoke_minimise<NUM>(_pFkt, vecParamNames, vecVals, vecErrs, vecFixed)
	bool bOk = 0;

	// stating all needed specialisations of the fit template function
	if(iNumParams == 1) bOk = __CALL_MINI(1);
	else if(iNumParams == 2) bOk = __CALL_MINI(2);
	else if(iNumParams == 3) bOk = __CALL_MINI(3);
	else if(iNumParams == 4) bOk = __CALL_MINI(4);
	else if(iNumParams == 5) bOk = __CALL_MINI(5);
	else if(iNumParams == 6) bOk = __CALL_MINI(6);
	else if(iNumParams == 7) bOk = __CALL_MINI(7);
	else if(iNumParams == 8) bOk = __CALL_MINI(8);
	else if(iNumParams == 9) bOk = __CALL_MINI(9);
	else if(iNumParams == 10) bOk = __CALL_MINI(10);
	else if(iNumParams == 11) bOk = __CALL_MINI(11);
	else if(iNumParams == 12) bOk = __CALL_MINI(12);
	else if(iNumParams == 13) bOk = __CALL_MINI(13);
	else if(iNumParams == 14) bOk = __CALL_MINI(14);
	else if(iNumParams == 15) bOk = __CALL_MINI(15);
	else tl::log_err("In ", __func__, ": Invalid number of function arguments.");


	// write back fitted values & errors
	for(std::size_t iParam=0; iParam<iNumParams; ++iParam)
	{
		pValues[iParam] = vecVals[iParam];
		pErrors[iParam] = vecErrs[iParam];
	}


	return int(bOk);
}
