/**
 * tlibs2
 * fitting and interpolation library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2012-2021
 * @note Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_FITTER_H__
#define __TLIBS2_FITTER_H__

#if __has_include(<Minuit2/FCNBase.h>) && __has_include(<Minuit2/MnTraceObject.h>)
	#include <Minuit2/FCNBase.h>
	#include <Minuit2/MnFcn.h>
	#include <Minuit2/FunctionMinimum.h>
	#include <Minuit2/MnMigrad.h>
	#include <Minuit2/MnPrint.h>

	#define __TLIBS2_USE_MINUIT__
#else
	#pragma message("tlibs2: Disabling Minuit library (not found).")
#endif

#include <vector>
#include <span>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <limits>
#include <type_traits>

#include "log.h"
#include "expr.h"

#include "maths.h"


namespace tl2 {


// ----------------------------------------------------------------------------
// Minuit interface
// @see http://seal.cern.ch/documents/minuit/mnusersguide.pdf
// ----------------------------------------------------------------------------
#ifdef __TLIBS2_USE_MINUIT__
using t_real_min = std::invoke_result_t<
	decltype(&ROOT::Minuit2::MnFcn::Up), ROOT::Minuit2::MnFcn>;


// ----------------------------------------------------------------------------
// function models

template<class t_real>
class FitterFuncModel
{
public:
	virtual ~FitterFuncModel() = default;

	virtual bool SetParams(const std::vector<t_real>& vecParams) = 0;
	virtual t_real operator()(t_real x) const = 0;
	virtual FitterFuncModel<t_real>* copy() const = 0;
};


/**
 * interface using supplied functions
 * iNumArgs also includes the "x" parameter to the function, m_vecVals does not
 */
template<class t_real, std::size_t iNumArgs, typename t_func>
class FitterLamFuncModel : public FitterFuncModel<t_real>
{
protected:
	t_func m_func{};
	std::vector<t_real> m_vecVals{};
	bool m_bSeparateFreeParam = 1;	// separate "x" from parameters (for fitter)

public:
	FitterLamFuncModel(t_func func, bool bSeparateX=1)
		: m_func{func}, m_vecVals{}, m_bSeparateFreeParam{bSeparateX}
	{
		m_vecVals.resize(m_bSeparateFreeParam ? iNumArgs-1 : iNumArgs);
	}


	virtual bool SetParams(const std::vector<t_real>& vecParams) override
	{
		for(std::size_t i=0; i<std::min(vecParams.size(), m_vecVals.size()); ++i)
			m_vecVals[i] = vecParams[i];
		return true;
	}


	virtual t_real operator()(t_real x = t_real(0)) const override
	{
		std::vector<t_real> vecValsWithX;
		if(m_bSeparateFreeParam)
		{
			vecValsWithX.push_back(x);
			for(t_real d : m_vecVals) vecValsWithX.push_back(d);
		}

		const std::vector<t_real> *pvecVals = m_bSeparateFreeParam ? &vecValsWithX : &m_vecVals;
		t_real funcval = call<iNumArgs, t_func, t_real, std::vector>(m_func, *pvecVals);
		return funcval;
	}


	virtual FitterLamFuncModel* copy() const override
	{
		FitterLamFuncModel<t_real, iNumArgs, t_func>* pMod =
			new FitterLamFuncModel<t_real, iNumArgs, t_func>(m_func);

		pMod->m_vecVals = this->m_vecVals;
		pMod->m_bSeparateFreeParam = this->m_bSeparateFreeParam;

		return pMod;
	}
};


/**
 * interface using supplied functions
 * iNumArgs also includes the "x" parameter to the function, m_vecVals does not
 */
template<class t_real>
class FitterParsedFuncModel : public FitterFuncModel<t_real>
{
protected:
	std::string m_func;

	std::string m_xName = "x";
	const std::vector<std::string>& m_vecNames;
	std::vector<t_real> m_vecVals;

	ExprParser<t_real> m_expr{};


public:
	FitterParsedFuncModel(const std::string& func, const std::string& xName, const std::vector<std::string>& vecNames)
		: m_func{func}, m_xName{xName}, m_vecNames{vecNames}
	{
		if(!m_expr.parse(m_func))
			throw std::runtime_error("Could not parse function.");
	}


	virtual bool SetParams(const std::vector<t_real>& vecParams) override
	{
		m_vecVals.resize(vecParams.size());
		for(std::size_t i=0; i<std::min(vecParams.size(), m_vecVals.size()); ++i)
			m_vecVals[i] = vecParams[i];
		return true;
	}


	virtual t_real operator()(t_real x = t_real(0)) const override
	{
		// copy the parsed expression to be thread safe
		ExprParser<t_real> expr = m_expr;

		// x is not used for minimiser
		if(m_xName != "")
			expr.register_var(m_xName, x);

		for(std::size_t i=0; i<m_vecVals.size(); ++i)
			expr.register_var(m_vecNames[i], m_vecVals[i]);

		t_real val = expr.eval();
		//std::cout << "f(" << x << ") = " << val << std::endl;
		return val;
	}


	virtual FitterParsedFuncModel* copy() const override
	{
		return new FitterParsedFuncModel<t_real>(m_func, m_xName, m_vecNames);
	}
};

// ----------------------------------------------------------------------------



/**
 * generic chi^2 calculation for fitting
 */
template<class t_real = t_real_min>
class Chi2Function : public ROOT::Minuit2::FCNBase
{
protected:
	const FitterFuncModel<t_real_min> *m_pfkt = nullptr;

	std::size_t m_num_pts = 0;
	const t_real* m_px = nullptr;
	const t_real* m_py = nullptr;
	const t_real* m_pdy = nullptr;

	t_real_min m_dSigma = 1.;
	bool m_bDebug = 0;


public:
	Chi2Function(const FitterFuncModel<t_real_min>* fkt=0,
		std::size_t num_pts=0, const t_real *px=0,
		const t_real *py=0, const t_real *pdy=0)
		: m_pfkt{fkt}, m_num_pts{num_pts}, m_px{px}, m_py{py}, m_pdy{pdy}
	{}

	virtual ~Chi2Function() = default;


	const Chi2Function<t_real> operator=(const Chi2Function<t_real>& other)
	{
		this->m_pfkt = other.m_pfkt;
		this->m_px = other.m_px;
		this->m_py = other.m_py;
		this->m_pdy = other.m_pdy;
		this->m_dSigma = other.m_dSigma;
		this->m_bDebug = other.m_bDebug;
	}

	Chi2Function(const Chi2Function<t_real>& other)
	{
		operator=(other);
	}


	/*
	 * chi^2 calculation
	 * based on the example in the Minuit user's guide:
	 * http://seal.cern.ch/documents/minuit/mnusersguide.pdf
	 */
	t_real_min chi2(const std::vector<t_real_min>& vecParams) const
	{
		// cannot operate on m_pfkt directly, because Minuit
		// uses more than one thread!
		std::unique_ptr<FitterFuncModel<t_real_min>> uptrFkt(m_pfkt->copy());
		FitterFuncModel<t_real_min>* pfkt = uptrFkt.get();

		pfkt->SetParams(vecParams);
		return tl2::chi2<t_real_min, decltype(*pfkt), const t_real*>(*pfkt, m_num_pts, m_px, m_py, m_pdy);
	}

	virtual t_real_min Up() const override { return m_dSigma*m_dSigma; }

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		t_real_min dChi2 = chi2(vecParams);
		if(m_bDebug) log_debug("chi2 = ", dChi2);
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
	const FitterFuncModel<t_real_min> *m_pfkt = nullptr;
	t_real_min m_dSigma = 1.;

public:
	MiniFunction(const FitterFuncModel<t_real_min>* fkt=0) : m_pfkt(fkt) {}
	virtual ~MiniFunction() = default;

	MiniFunction(const MiniFunction<t_real>& other)
		: m_pfkt(other.m_pfkt), m_dSigma(other.m_dSigma)
	{}

	const MiniFunction<t_real> operator=(const MiniFunction<t_real>& other)
	{
		this->m_pfkt = other.m_pfkt;
		this->m_dSigma = other.m_dSigma;
	}

	virtual t_real_min Up() const override { return m_dSigma*m_dSigma; }

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		// cannot operate on m_pfkt directly, because Minuit
		// uses more than one thread!
		std::unique_ptr<FitterFuncModel<t_real_min>> uptrFkt(m_pfkt->copy());
		FitterFuncModel<t_real_min>* pfkt = uptrFkt.get();

		pfkt->SetParams(vecParams);
		return (*pfkt)(t_real_min(0));	// "0" is an ignored dummy value here
	}

	void SetSigma(t_real_min dSig) { m_dSigma = dSig; }
	t_real_min GetSigma() const { return m_dSigma; }
};



// ----------------------------------------------------------------------------


/**
 * fit function to x,y,dy data points
 */
template<class t_real = t_real_min, std::size_t iNumArgs, typename t_func>
bool fit(t_func&& func,

	const std::vector<t_real>& vecX,
	const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr,

	const std::vector<std::string>& vecParamNames,	// size: iNumArgs-1
	std::vector<t_real>& vecVals,
	std::vector<t_real>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr,

	bool bDebug=1) noexcept
{
	try
	{
		if(!vecX.size() || !vecY.size() || !vecYErr.size())
		{
			log_err("No data given to fitter.");
			return false;
		}

		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b)->bool { return b; }))
			{
				log_err("All parameters are fixed.");
				return false;
			}

		// convert vectors if value types don't match with minuit's type
		std::vector<t_real_min> vecXConverted, vecYConverted, vecYErrConverted;
		if constexpr(!std::is_same_v<t_real, t_real_min>)
		{
			vecXConverted.reserve(vecX.size());
			vecYConverted.reserve(vecY.size());
			vecYErrConverted.reserve(vecYErr.size());

			for(t_real d : vecX) vecXConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecY) vecYConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecYErr) vecYErrConverted.push_back(static_cast<t_real_min>(d));
		}

		FitterLamFuncModel<t_real_min, iNumArgs, t_func> mod(func);

		std::unique_ptr<Chi2Function<t_real_min>> chi2;
		if constexpr(std::is_same_v<t_real, t_real_min>)
			chi2 = std::make_unique<Chi2Function<t_real_min>>(&mod, vecX.size(), vecX.data(), vecY.data(), vecYErr.data());
		else if constexpr(!std::is_same_v<t_real, t_real_min>)
			chi2 = std::make_unique<Chi2Function<t_real_min>>(&mod, vecXConverted.size(), vecXConverted.data(), vecYConverted.data(), vecYErrConverted.data());

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			params.Add(vecParamNames[iParam], static_cast<t_real_min>(vecVals[iParam]), static_cast<t_real_min>(vecErrs[iParam]));
			if(pVecFixed && (*pVecFixed)[iParam])
				params.Fix(vecParamNames[iParam]);
		}

		ROOT::Minuit2::MnMigrad migrad(*chi2, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bValidFit = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			vecVals[iParam] = static_cast<t_real>(mini.UserState().Value(vecParamNames[iParam]));
			vecErrs[iParam] = static_cast<t_real>(std::fabs(mini.UserState().Error(vecParamNames[iParam])));
		}

		if(bDebug)
			log_debug(mini);

		return bValidFit;
	}
	catch(const std::exception& ex)
	{
		log_err(ex.what());
	}

	return false;
}


/**
 * fit expression to x,y,dy data points
 */
template<class t_real = t_real_min>
bool fit_expr(const std::string& func,

	const std::vector<t_real>& vecX,
	const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr,

	const std::string& strXName,
	const std::vector<std::string>& vecParamNames,	// size: iNumArgs-1
	std::vector<t_real>& vecVals,
	std::vector<t_real>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr,

	bool bDebug=1) noexcept
{
	try
	{
		if(!vecX.size() || !vecY.size() || !vecYErr.size())
		{
			log_err("No data given to fitter.");
			return false;
		}

		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b)->bool { return b; }))
			{
				log_err("All parameters are fixed.");
				return false;
			}

		// convert vectors if value types don't match with minuit's type
		std::vector<t_real_min> vecXConverted, vecYConverted, vecYErrConverted;
		if constexpr(!std::is_same_v<t_real, t_real_min>)
		{
			vecXConverted.reserve(vecX.size());
			vecYConverted.reserve(vecY.size());
			vecYErrConverted.reserve(vecYErr.size());

			for(t_real d : vecX) vecXConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecY) vecYConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecYErr) vecYErrConverted.push_back(static_cast<t_real_min>(d));
		}

		FitterParsedFuncModel<t_real_min> mod(func, strXName, vecParamNames);

		std::unique_ptr<Chi2Function<t_real_min>> chi2;
		if constexpr(std::is_same_v<t_real, t_real_min>)
			chi2 = std::make_unique<Chi2Function<t_real_min>>(&mod, vecX.size(), vecX.data(), vecY.data(), vecYErr.data());
		else if constexpr(!std::is_same_v<t_real, t_real_min>)
			chi2 = std::make_unique<Chi2Function<t_real_min>>(&mod, vecXConverted.size(), vecXConverted.data(), vecYConverted.data(), vecYErrConverted.data());

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			params.Add(vecParamNames[iParam], static_cast<t_real_min>(vecVals[iParam]), static_cast<t_real_min>(vecErrs[iParam]));
			if(pVecFixed && (*pVecFixed)[iParam])
				params.Fix(vecParamNames[iParam]);
		}

		ROOT::Minuit2::MnMigrad migrad(*chi2, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bValidFit = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			vecVals[iParam] = static_cast<t_real>(mini.UserState().Value(vecParamNames[iParam]));
			vecErrs[iParam] = static_cast<t_real>(std::fabs(mini.UserState().Error(vecParamNames[iParam])));
		}

		if(bDebug)
			log_debug(mini);

		return bValidFit;
	}
	catch(const std::exception& ex)
	{
		log_err(ex.what());
	}

	return false;
}


/**
 * find function minimum
 */
template<class t_real=t_real_min, std::size_t iNumArgs, typename t_func>
bool minimise(t_func&& func, const std::vector<std::string>& vecParamNames,
	std::vector<t_real>& vecVals, std::vector<t_real>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr, bool bDebug=1) noexcept
{
	try
	{
		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b)->bool { return b; }))
			{
				log_err("All parameters are fixed.");
				return false;
			}

		FitterLamFuncModel<t_real_min, iNumArgs, t_func> mod(func, false);
		MiniFunction<t_real_min> chi2(&mod);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			params.Add(vecParamNames[iParam], static_cast<t_real_min>(vecVals[iParam]), static_cast<t_real_min>(vecErrs[iParam]));
			if(pVecFixed && (*pVecFixed)[iParam])
				params.Fix(vecParamNames[iParam]);
		}

		ROOT::Minuit2::MnMigrad migrad(chi2, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bMinimumValid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			vecVals[iParam] = static_cast<t_real>(mini.UserState().Value(vecParamNames[iParam]));
			vecErrs[iParam] = static_cast<t_real>(std::fabs(mini.UserState().Error(vecParamNames[iParam])));
		}

		if(bDebug)
			log_debug(mini);

		return bMinimumValid;
	}
	catch(const std::exception& ex)
	{
		log_err(ex.what());
	}

	return false;
}


/**
 * find function minimum for an expression
 */
template<class t_real = t_real_min>
bool minimise_expr(const std::string& func, const std::vector<std::string>& vecParamNames,
	std::vector<t_real>& vecVals, std::vector<t_real>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr, bool bDebug=1) noexcept
{
	try
	{
		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b)->bool { return b; }))
			{
				log_err("All parameters are fixed.");
				return false;
			}

		FitterParsedFuncModel<t_real_min> mod(func, "", vecParamNames);
		MiniFunction<t_real_min> chi2(&mod);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			params.Add(vecParamNames[iParam], static_cast<t_real_min>(vecVals[iParam]), static_cast<t_real_min>(vecErrs[iParam]));
			if(pVecFixed && (*pVecFixed)[iParam])
				params.Fix(vecParamNames[iParam]);
		}

		ROOT::Minuit2::MnMigrad migrad(chi2, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bMinimumValid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			vecVals[iParam] = static_cast<t_real>(mini.UserState().Value(vecParamNames[iParam]));
			vecErrs[iParam] = static_cast<t_real>(std::fabs(mini.UserState().Error(vecParamNames[iParam])));
		}

		if(bDebug)
			log_debug(mini);

		return bMinimumValid;
	}
	catch(const std::exception& ex)
	{
		log_err(ex.what());
	}

	return false;
}

#endif	// __TLIBS2_USE_MINUIT__
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// interpolation
// ----------------------------------------------------------------------------

/**
 * @see http://mathworld.wolfram.com/BernsteinPolynomial.html
 */
template<typename T> T bernstein(int i, int n, T t)
{
	T bino = boost::math::binomial_coefficient<T>(n, i);
	return bino * pow(t, i) * pow(1-t, n-i);
}

/**
 * @see http://mathworld.wolfram.com/BezierCurve.html
 */
template<class t_vec, typename T=typename t_vec::value_type>
t_vec bezier(const t_vec* P, std::size_t N, T t)
{
	if(N==0) return t_vec{};
	const int n = N-1;

	t_vec vec(P[0].size());
	for(std::size_t i=0; i<vec.size(); ++i) vec[i] = T(0);

	for(int i=0; i<=n; ++i)
		vec += P[i]*bernstein(i, n, t);

	return vec;
}


/**
 * @see http://mathworld.wolfram.com/B-Spline.html
 */
template<typename T>
T bspline_base(int i, int j, T t, const std::vector<T>& knots)
{
	if(j==0)
	{
		if((knots[i] <= t) && (t < knots[i+1]) && (knots[i]<knots[i+1]))
			return 1.;
		return 0.;
	}

	T val11 = (t - knots[i]) / (knots[i+j]-knots[i]);
	T val12 = bspline_base(i, j-1, t, knots);
	T val1 = val11 * val12;

	T val21 = (knots[i+j+1]-t) / (knots[i+j+1]-knots[i+1]);
	T val22 = bspline_base(i+1, j-1, t, knots);
	T val2 = val21 * val22;

	T val = val1 + val2;
	return val;
}


/**
 * @see http://mathworld.wolfram.com/B-Spline.html
 */
template<class t_vec, typename T=typename t_vec::value_type>
t_vec bspline(const t_vec* P, std::size_t N, T t, const std::vector<T>& knots)
{
	if(N==0) return t_vec{};
	const int n = N-1;
	const int m = knots.size()-1;
	const int degree = m-n-1;

	t_vec vec(P[0].size());
	for(std::size_t i=0; i<vec.size(); ++i)
		vec[i] = T(0);

	for(int i=0; i<=n; ++i)
		vec += P[i]*bspline_base(i, degree, t, knots);

	return vec;
}


// ----------------------------------------------------------------------------

template<class t_vec, typename T=typename t_vec::value_type>
class Bezier
{
	protected:
		std::unique_ptr<t_vec[]> m_pvecs;
		std::size_t m_iN;


	public:
		Bezier(std::size_t N, const T *px, const T *py) : m_iN(N)
		{
			m_pvecs.reset(new t_vec[m_iN]);

			for(std::size_t i=0; i<m_iN; ++i)
			{
				m_pvecs[i].resize(2);
				m_pvecs[i][0] = px[i];
				m_pvecs[i][1] = py[i];
			}
		}


		t_vec operator()(T t) const
		{
			return bezier<t_vec, T>(m_pvecs.get(), m_iN, t);
		}
};


template<class t_vec, typename T=typename t_vec::value_type>
class BSpline
{
	protected:
		T m_eps = std::numeric_limits<T>::epsilon();
		std::unique_ptr<t_vec[]> m_pvecs{};
		std::size_t m_iN, m_iDegree{};
		std::vector<T> m_vecKnots{};


	public:
		BSpline(std::size_t N, const T *px, const T *py,
			unsigned int iDegree=3) : m_iN(N), m_iDegree(iDegree)
		{
			m_pvecs.reset(new t_vec[m_iN]);

			for(std::size_t i=0; i<m_iN; ++i)
			{
				m_pvecs[i].resize(2);
				m_pvecs[i][0] = px[i];
				m_pvecs[i][1] = py[i];
			}

			std::size_t iM = m_iDegree + m_iN + 1;
			m_vecKnots.resize(iM);


			// set knots to uniform, nonperiodic B-Spline
			for(unsigned int i=0; i<m_iDegree+1; ++i)
				m_vecKnots[i] = 0.+i*m_eps;
			for(unsigned int i=iM-m_iDegree-1; i<iM; ++i)
				m_vecKnots[i] = 1.-i*m_eps;
			for(unsigned int i=m_iDegree+1; i<iM-m_iDegree-1; ++i)
				m_vecKnots[i] = T(i+1-m_iDegree-1) / T(iM-2*m_iDegree-2 + 1);
		}


		t_vec operator()(T t) const
		{
			if(m_iN==0)
			{
				t_vec vecNull(2);
				vecNull[0] = vecNull[1] = 0.;
				return vecNull;
			}

			t_vec vec = bspline<t_vec, T>(m_pvecs.get(), m_iN, t, m_vecKnots);

			// remove epsilon dependence
			if(t<=0.) vec = m_pvecs[0];
			if(t>=1.) vec = m_pvecs[m_iN-1];

			return vec;
		}


		void SetEps(T eps) { m_eps = eps; }
};


template<class t_vec, typename T=typename t_vec::value_type>
class LinInterp
{
protected:
	std::unique_ptr<t_vec[]> m_pvecs;
	std::size_t m_iN;

public:
	LinInterp(std::size_t N, const T *px, const T *py) : m_iN(N)
	{
		m_pvecs.reset(new t_vec[m_iN]);

		for(std::size_t i=0; i<m_iN; ++i)
		{
			m_pvecs[i].resize(2);
			m_pvecs[i][0] = px[i];
			m_pvecs[i][1] = py[i];
		}

		// ensure that vector is sorted by x values
		std::stable_sort(m_pvecs.get(), m_pvecs.get()+m_iN,
			[](const t_vec& vec1, const t_vec& vec2) -> bool
			{ return vec1[0] < vec2[0]; });
	}


	T operator()(T x) const
	{
		const auto* iterBegin = m_pvecs.get();
		const auto* iterEnd = m_pvecs.get() + m_iN;

		if(m_iN == 0) return T(0);
		if(m_iN == 1) return (*iterBegin)[1];

		const auto* iterLower = std::lower_bound(iterBegin, iterEnd, x,
			[](const t_vec& vec, const T& x) -> bool
			{ return vec[0] < x; });

		// lower bound at end of range?
		if(iterLower == iterEnd || iterLower == iterEnd-1)
			iterLower = iterEnd - 2;
		const auto* iter2 = iterLower + 1;

		T xrange = (*iter2)[0] - (*iterLower)[0];
		T xpos = (x-(*iterLower)[0]) / xrange;

		return std::lerp((*iterLower)[1], (*iter2)[1], xpos);
	}
};
// ------------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------------
// peak finder
// ------------------------------------------------------------------------------------------------
template<class T, std::size_t N> struct __sort_obj
{
	T vec[N];
};


template<class T, std::size_t N>
bool __comp_fkt(const __sort_obj<T, N>& t0, const __sort_obj<T, N>& t1)
{
	return t0.vec[0] < t1.vec[0];
}


/**
 * simultaneously sort two arrays
 */
template<class Iter = double*>
void __sort_2(Iter begin1, Iter end1, Iter begin2)
{
	using T = typename std::iterator_traits<Iter>::value_type;

	const std::size_t N = end1 - begin1;
	std::unique_ptr<__sort_obj<T, 2>[]>obj{new __sort_obj<T, 2>[N]};

	for(std::size_t i=0; i<N; ++i)
	{
		obj.get()[i].vec[0] = *(begin1+i);
		obj.get()[i].vec[1] = *(begin2+i);
	}

	std::stable_sort(obj.get(), obj.get()+N, __comp_fkt<T, 2>);

	for(std::size_t i=0; i<N; ++i)
	{
		*(begin1+i) = obj.get()[i].vec[0];
		*(begin2+i) = obj.get()[i].vec[1];
	}
}


/**
 * find the zeros of a curve
 */
template<class t_cont>
std::vector<std::size_t> find_zeroes(const t_cont& cont)
{
	const std::size_t num_pts = cont.size();
	std::vector<std::size_t> indices;

	for(std::size_t i=0; i<num_pts-1; ++i)
	{
		if(std::signbit(cont[i]) != std::signbit(cont[i+1]))
			indices.push_back(i);
	}

	return indices;
}


/**
 * try to identify local minima and maxima on a curve
 */
template<typename T>
void find_peaks(std::size_t num_pts, const T* _px, const T* _py, unsigned int spline_order,
	std::vector<T>& peaks_x, std::vector<T>& peaks_sizes,
	std::vector<T>& peaks_widths, std::vector<bool>& peaks_minima,
	std::size_t num_spline = 512, T eps = std::numeric_limits<T>::epsilon())
{
	if(num_pts < 1)
		return;

	using t_vec = tl2::vec<T, std::vector>;

	// allocate memory
	std::unique_ptr<T[]> uptrMem{new T[2*num_pts + 2*num_spline]};
	T *px = uptrMem.get() + 0*num_pts;
	T *py = uptrMem.get() + 1*num_pts;
	T *spline_x = uptrMem.get() + 2*num_pts + 0*num_spline;
	T *spline_y = uptrMem.get() + 2*num_pts + 1*num_spline;


	// sort input values by x
	std::copy(_px, _px+num_pts, px);
	std::copy(_py, _py+num_pts, py);
	__sort_2<T*>(px, px+num_pts, py);


	BSpline<t_vec, T> spline(num_pts, px, py, spline_order);
	spline.SetEps(eps);
	const T* y_min = std::min_element(py, py+num_pts);

	for(std::size_t iSpline=0; iSpline<num_spline; ++iSpline)
	{
		const T dT = T(iSpline) / T(num_spline-1);
		t_vec vec = spline(dT);

		spline_x[iSpline] = vec[0];
		spline_y[iSpline] = vec[1];
	}


	// differentiate curve
	std::vector<T> spline_diff = tl2::diff<std::vector<T>>(
		std::span<T>(spline_x, num_spline), std::span<T>(spline_y, num_spline));
	std::vector<T> spline_diff2 = tl2::diff<std::vector<T>>(
		std::span<T>(spline_x, num_spline), std::span<T>(spline_diff.data(), num_spline));

	// find zeros
	std::vector<std::size_t> zeros = find_zeroes(spline_diff);


	for(std::size_t zero_idx=0; zero_idx<zeros.size(); ++zero_idx)
	{
		const std::size_t cur_zero_idx = zeros[zero_idx];

		// saddle point
		if(tl2::equals_0<T>(spline_diff2[cur_zero_idx], eps))
			continue;

		// minima or maxima
		bool is_minimum = (spline_diff2[cur_zero_idx] > 0.);

		int min_idx_left = -1;
		int min_idx_right = -1;
		if(zero_idx > 0)
			min_idx_left = zeros[zero_idx-1];
		if(zero_idx < zeros.size() - 1)
			min_idx_right = zeros[zero_idx+1];

		T height = 0.;
		T width = 0.;
		T div = 0.;

		// extremum left of the peak
		if(min_idx_left >= 0)
		{
			height += (spline_y[cur_zero_idx] - spline_y[min_idx_left]);
			width += std::abs((spline_x[cur_zero_idx] - spline_x[min_idx_left]));
			div += 1.;
		}

		// extremum right of the peak
		if(min_idx_right >= 0)
		{
			height += (spline_y[cur_zero_idx] - spline_y[min_idx_right]);
			width += std::abs((spline_x[cur_zero_idx] - spline_x[min_idx_right]));
			div += 1.;
		}

		// no adjacent extrema...
		if(min_idx_left < 0 && min_idx_right < 0)
		{
			height = spline_y[cur_zero_idx] - *y_min;
			width = (px[num_pts-1] - px[0]) / 10.;	// guess something...
			div = 1.;
		}

		if(div != 0.)
		{
			height /= div;
			width /= div;
		}

		peaks_x.push_back(spline_x[cur_zero_idx]);
		peaks_sizes.push_back(height);
		peaks_widths.push_back(width);
		peaks_minima.push_back(is_minimum);
	}
}
// ----------------------------------------------------------------------------


}

#endif
