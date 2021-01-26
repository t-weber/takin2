/**
 * Minuit interface
 *
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date April 2012
 * @license GPLv2 or GPLv3
 *
 * @desc general fitter structure (i.e. function => chi^2 calculation => calling
 * 	minuit) originally based on the examples in the Minuit user's guide:
 * 	http://seal.cern.ch/documents/minuit/mnusersguide.pdf
 */

#ifndef __MINUIT_IFACE_H__
#define __MINUIT_IFACE_H__

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnFcn.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <type_traits>
#include <limits>
#include <cmath>

#include "../helper/boost_hacks.h"
#include <boost/numeric/ublas/vector.hpp>

#include "funcmod.h"
#include "../helper/misc.h"
#include "../log/log.h"
#include "../math/stat.h"


namespace tl {
//using t_real_min = double;
using t_real_min = typename std::result_of<
	decltype(&ROOT::Minuit2::MnFcn::Up)(ROOT::Minuit2::MnFcn)>::type;

// function models
using MinuitFuncModel = class FitterFuncModel<t_real_min>;
using MinuitFuncModel_nd = class FitterFuncModel_nd<t_real_min>;

template<class t_real>
using MinuitMultiFuncModel = class FitterMultiFuncModel<t_real, t_real_min>;

template<std::size_t iNumArgs, typename t_func>
using MinuitLamFuncModel = FitterLamFuncModel<t_real_min, iNumArgs, t_func>;


/**
 * generic chi^2 calculation for fitting
 */
template<class t_real = t_real_min>
class Chi2Function : public ROOT::Minuit2::FCNBase
{
protected:
	const MinuitFuncModel *m_pfkt = nullptr;

	std::size_t m_uiLen;
	const t_real* m_px = nullptr;
	const t_real* m_py = nullptr;
	const t_real* m_pdy = nullptr;

	t_real_min m_dSigma = 1.;
	bool m_bDebug = 0;

public:
	Chi2Function(const MinuitFuncModel* fkt=0,
		std::size_t uiLen=0, const t_real *px=0,
		const t_real *py=0, const t_real *pdy=0)
		: m_pfkt(fkt), m_uiLen(uiLen), m_px(px), m_py(py), m_pdy(pdy)
	{}

	virtual ~Chi2Function() = default;

	/*
	 * chi^2 calculation
	 * based on the example in the Minuit user's guide:
	 * http://seal.cern.ch/documents/minuit/mnusersguide.pdf
	 */
	t_real_min chi2(const std::vector<t_real_min>& vecParams) const
	{
		// cannot operate on m_pfkt directly, because Minuit
		// uses more than one thread!
		std::unique_ptr<MinuitFuncModel> uptrFkt(m_pfkt->copy());
		MinuitFuncModel* pfkt = uptrFkt.get();

		/*bool bParamsOk = */pfkt->SetParams(vecParams);
		//if(!bParamsOk)
		//	return std::numeric_limits<t_real_min>::max();

		return tl::chi2<t_real_min, decltype(*pfkt), const t_real*>(*pfkt, m_uiLen, m_px, m_py, m_pdy);
	}

	virtual t_real_min Up() const override { return m_dSigma*m_dSigma; }

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		t_real_min dChi2 = chi2(vecParams);
		if(m_bDebug) tl::log_debug("Chi2 = ", dChi2);
		return dChi2;
	}

	void SetSigma(t_real_min dSig) { m_dSigma = dSig; }
	t_real_min GetSigma() const { return m_dSigma; }

	void SetDebug(bool b) { m_bDebug = b; }
};


/**
 * function adaptor for minimisation
 */
template<class t_real = t_real_min>
class MiniFunction : public ROOT::Minuit2::FCNBase
{
protected:
	const MinuitFuncModel *m_pfkt = nullptr;
	t_real_min m_dSigma = 1.;

public:
	MiniFunction(const MinuitFuncModel* fkt=0) : m_pfkt(fkt) {}
	virtual ~MiniFunction() = default;

	virtual t_real_min Up() const override { return m_dSigma*m_dSigma; }

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		// cannot operate on m_pfkt directly, because Minuit
		// uses more than one thread!
		std::unique_ptr<MinuitFuncModel> uptrFkt(m_pfkt->copy());
		MinuitFuncModel* pfkt = uptrFkt.get();

		pfkt->SetParams(vecParams);
		return (*pfkt)(t_real_min(0));	// "0" is an ignored dummy value here
	}

	void SetSigma(t_real_min dSig) { m_dSigma = dSig; }
	t_real_min GetSigma() const { return m_dSigma; }
};


// ----------------------------------------------------------------------------


/**
 * chi^2 calculation using multiple simultaneous functions where each
 * function can additionally have multiple parameter sets
 */
template<class t_real = t_real_min, template<class...> class t_cont=std::vector>
class Chi2Function_mult : public ROOT::Minuit2::FCNBase
{
protected:
	t_cont<const MinuitMultiFuncModel<t_real>*> m_vecFkt;

	t_cont<std::size_t> m_vecLen;
	t_cont<const t_real*> m_vecpX;
	t_cont<const t_real*> m_vecpY;
	t_cont<const t_real*> m_vecpDY;

	t_real_min m_dSigma = 1.;
	bool m_bDebug = 0;

public:
	Chi2Function_mult() = default;
	virtual ~Chi2Function_mult() = default;

	void AddFunc(const MinuitMultiFuncModel<t_real>* pMod, std::size_t iNumDat,
		const t_real *pX, const t_real *pY, const t_real *pdY)
	{
		m_vecFkt.push_back(pMod);
		m_vecLen.push_back(iNumDat);
		m_vecpX.push_back(pX);
		m_vecpY.push_back(pY);
		m_vecpDY.push_back(pdY);
	}

	t_real_min chi2(std::size_t iFkt, const std::vector<t_real_min>& vecParams) const
	{
		std::unique_ptr<MinuitMultiFuncModel<t_real>> uptrFkt(m_vecFkt[iFkt]->copy());
		MinuitMultiFuncModel<t_real>* pfkt = uptrFkt.get();

		const std::size_t iNumParamSets = pfkt->GetParamSetCount();
		t_real_min dChi = t_real_min(0);
		for(std::size_t iParamSet=0; iParamSet<iNumParamSets; ++iParamSet)
		{
			pfkt->SetParamSet(iParamSet);
			std::size_t iLen = pfkt->GetExpLen();
			const t_real* pX = pfkt->GetExpX();
			const t_real* pY = pfkt->GetExpY();
			const t_real* pDY = pfkt->GetExpDY();

			// default experimental values if none are given in the model
			if(!pX || !pY || !pDY)
			{
				iLen = m_vecLen[iFkt];
				pX = m_vecpX[iFkt];
				pY = m_vecpY[iFkt];
				pDY = m_vecpDY[iFkt];
			}

			pfkt->SetParams(vecParams);

			t_real_min dSingleChi = tl::chi2<t_real_min, decltype(*pfkt), const t_real*>
				(*pfkt, iLen, pX, pY, pDY);
			dChi += dSingleChi;

			if(m_bDebug && iNumParamSets>1)
				tl::log_debug("Function ", iParamSet, " chi2 = ", dSingleChi);
		}
		dChi /= t_real_min(iNumParamSets);
		return dChi;
	}

	t_real_min chi2(const std::vector<t_real_min>& vecParams) const
	{
		t_real_min dChi = t_real_min(0);
		for(std::size_t iFkt=0; iFkt<m_vecFkt.size(); ++iFkt)
		{
			const t_real_min dSingleChi = chi2(iFkt, vecParams);
			dChi += dSingleChi;

			if(m_bDebug && m_vecFkt.size()>1)
				tl::log_debug("Function ", iFkt, " chi2 = ", dSingleChi);
		}
		dChi /= t_real_min(m_vecFkt.size());
		return dChi;
	}

	virtual t_real_min Up() const override { return m_dSigma*m_dSigma; }

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		t_real_min dChi2 = chi2(vecParams);
		if(m_bDebug) tl::log_debug("Total chi2 = ", dChi2);
		return dChi2;
	}

	void SetSigma(t_real_min dSig) { m_dSigma = dSig; }
	t_real_min GetSigma() const { return m_dSigma; }

	void SetDebug(bool b) { m_bDebug = b; }
};

// for most cases data type of measured values and internal data type is the same: t_real_min
using Chi2FunctionMult = Chi2Function_mult<t_real_min, std::vector>;


// ----------------------------------------------------------------------------


// in n dimensions
template<class t_real = t_real_min, template<class...> class t_cont = std::vector>
class Chi2Function_nd : public ROOT::Minuit2::FCNBase
{
protected:
	const MinuitFuncModel_nd *m_pfkt = nullptr;

	const t_cont<t_cont<t_real>> *m_pvecvecX = nullptr;
	const t_cont<t_real> *m_pvecY = nullptr;
	const t_cont<t_real> *m_pvecDY = nullptr;

	bool m_bDebug = 0;

public:
	Chi2Function_nd(const MinuitFuncModel_nd* fkt,
		const t_cont<t_cont<t_real>> &vecvecX,
		const t_cont<t_real> &vecY, const t_cont<t_real> &vecDY)
		: m_pfkt(fkt), m_pvecvecX(&vecvecX), m_pvecY(&vecY), m_pvecDY(&vecDY)
	{}

	virtual ~Chi2Function_nd() = default;


	t_real_min chi2(const std::vector<t_real_min>& vecParams) const
	{
		std::unique_ptr<MinuitFuncModel_nd> uptrFkt(m_pfkt->copy());
		MinuitFuncModel_nd* pfkt = uptrFkt.get();
		pfkt->SetParams(vecParams);

		return tl::chi2_nd<t_real_min, t_real, decltype(*pfkt), t_cont>
			(*pfkt, *m_pvecvecX, *m_pvecY, *m_pvecDY);
	}

	virtual t_real_min Up() const override
	{
		// 1. for chi^2
		return 1.;
	}

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		t_real_min dChi2 = chi2(vecParams);
		if(m_bDebug) tl::log_debug("Chi2 = ", dChi2);
		return dChi2;
	}

	void SetDebug(bool b) { m_bDebug = b; }
};


// ----------------------------------------------------------------------------


/**
 * fit function to x,y,dy data points
 */
template<std::size_t iNumArgs, typename t_func>
bool fit(t_func&& func,

	const std::vector<t_real_min>& vecX,
	const std::vector<t_real_min>& vecY,
	const std::vector<t_real_min>& vecYErr,

	const std::vector<std::string>& vecParamNames,	// size: iNumArgs-1
	std::vector<t_real_min>& vecVals,
	std::vector<t_real_min>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr,

	bool bDebug=1) noexcept
{
	try
	{
		if(!vecX.size() || !vecY.size() || !vecYErr.size())
		{
			tl::log_err("No data given to fitter.");
			return false;
		}

		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b)->bool { return b; }))
			{
				tl::log_err("All parameters are fixed.");
				return false;
			}

		MinuitLamFuncModel<iNumArgs, t_func> mod(func);
		Chi2Function<t_real_min> chi2(&mod, vecX.size(), vecX.data(), vecY.data(), vecYErr.data());

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			params.Add(vecParamNames[iParam], vecVals[iParam], vecErrs[iParam]);
			if(pVecFixed && (*pVecFixed)[iParam])
				params.Fix(vecParamNames[iParam]);
		}

		ROOT::Minuit2::MnMigrad migrad(chi2, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bValidFit = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			vecVals[iParam] = mini.UserState().Value(vecParamNames[iParam]);
			vecErrs[iParam] = std::fabs(mini.UserState().Error(vecParamNames[iParam]));
		}

		if(bDebug)
			tl::log_debug(mini);

		return bValidFit;
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}

	return false;
}


/**
 * find function minimum
 */
template<std::size_t iNumArgs, typename t_func>
bool minimise(t_func&& func, const std::vector<std::string>& vecParamNames,
	std::vector<t_real_min>& vecVals, std::vector<t_real_min>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr, bool bDebug=1) noexcept
{
	try
	{
		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b)->bool { return b; }))
			{
				tl::log_err("All parameters are fixed.");
				return false;
			}

		MinuitLamFuncModel<iNumArgs, t_func> mod(func, false);
		MiniFunction<t_real_min> chi2(&mod);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			params.Add(vecParamNames[iParam], vecVals[iParam], vecErrs[iParam]);
			if(pVecFixed && (*pVecFixed)[iParam])
				params.Fix(vecParamNames[iParam]);
		}

		ROOT::Minuit2::MnMigrad migrad(chi2, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bMinimumValid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			vecVals[iParam] = mini.UserState().Value(vecParamNames[iParam]);
			vecErrs[iParam] = std::fabs(mini.UserState().Error(vecParamNames[iParam]));
		}

		if(bDebug)
			tl::log_debug(mini);

		return bMinimumValid;
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}

	return false;
}

// ----------------------------------------------------------------------------

}

#endif
