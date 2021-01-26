/**
 * abstract function model base classes
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013
 * @license GPLv2 or GPLv3
 */

#ifndef __FUNC_MOD_H__
#define __FUNC_MOD_H__

#include <string>
#include <vector>
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include "../helper/traits.h"

namespace tl {

/**
 * parametric function
 */
template<class t_vec = boost::numeric::ublas::vector<double>, class T = typename t_vec::value_type>
class FunctionModel_param
{
public:
	virtual ~FunctionModel_param() = default;

	// t = 0..1
	virtual t_vec operator()(T t) const = 0;
	virtual const char* GetModelName() const = 0;
};



// ----------------------------------------------------------------------------


/**
 * explicit function
 */
template<class T = double> class FunctionModel
{
public:
	virtual ~FunctionModel() = default;

	virtual T operator()(T x) const = 0;
	virtual const char* GetModelName() const = 0;
};

// synonym
//template<class T=double> using FunctionModel_gen = class FunctionModel<T>;


template<class t_real>
class FitterFuncModel : public FunctionModel<t_real>
{
public:
	virtual ~FitterFuncModel() = default;

	virtual bool SetParams(const std::vector<t_real>& vecParams) = 0;
	virtual t_real operator()(t_real x) const override = 0;

	virtual FitterFuncModel* copy() const = 0;
	virtual std::string print(bool bFillInSyms=true) const /*= 0;*/ { return ""; }

	virtual const char* GetModelName() const override /*= 0;*/ { return "FitterFuncModel"; }
	virtual std::vector<std::string> GetParamNames() const = 0;
	virtual std::vector<t_real> GetParamValues() const = 0;
	virtual std::vector<t_real> GetParamErrors() const = 0;

	friend std::ostream& operator<<(std::ostream& ostr, const FitterFuncModel<t_real>& fkt)
	{
		ostr << fkt.print();
		return ostr;
	}
};


// ----------------------------------------------------------------------------


/**
 * explicit function with multiple internal parameter sets
 */
template<class T = double> class FunctionModel_multi : public FunctionModel<T>
{
public:
	virtual ~FunctionModel_multi() = default;

	virtual std::size_t GetParamSetCount() const = 0;
	virtual void SetParamSet(std::size_t iSet) = 0;
};

// synonym
//template<class T=double> using FunctionModel_multi_gen = class FunctionModel_multi<T>;



template<class t_real, class t_real_min>
class FitterMultiFuncModel : public FunctionModel_multi<t_real_min>
{
public:
	virtual ~FitterMultiFuncModel() = default;

	virtual bool SetParams(const std::vector<t_real_min>& vecParams) = 0;
	virtual t_real_min operator()(t_real_min x) const override = 0;

	virtual FitterMultiFuncModel* copy() const = 0;
	virtual std::string print(bool bFillInSyms=true) const /*= 0;*/ { return ""; }

	virtual const char* GetModelName() const override /*= 0;*/ { return "FitterMultiFuncModel"; }
	virtual std::vector<std::string> GetParamNames() const = 0;
	virtual std::vector<t_real_min> GetParamValues() const = 0;
	virtual std::vector<t_real_min> GetParamErrors() const = 0;

	virtual void SetParamSet(std::size_t iSet) override /*= 0;*/ {}
	virtual std::size_t GetParamSetCount() const override /*= 0;*/ { return 1; }
	// optional intrinsic measured values for multi-parameter functions
	virtual std::size_t GetExpLen() const /*= 0;*/ { return 0; }
	virtual const t_real* GetExpX() const /*= 0;*/ { return nullptr; }
	virtual const t_real* GetExpY() const /*= 0*/ { return nullptr; }
	virtual const t_real* GetExpDY() const /*= 0*/ { return nullptr; }

	friend std::ostream& operator<<(std::ostream& ostr, const FitterMultiFuncModel<t_real, t_real_min>& fkt)
	{
		ostr << fkt.print();
		return ostr;
	}
};


// ----------------------------------------------------------------------------


/**
 * interface for n dimensional function
 */
template<class T = double> class FunctionModel_nd
{
protected:

public:
	virtual ~FunctionModel_nd() = default;

	virtual T operator()(const std::vector<T>& vecx) const = 0;
	virtual const char* GetModelName() const = 0;
};

// synonym
//template<class T=double> using FunctionModel_nd_gen = class FunctionModel_nd<T>;


template<class t_real>
class FitterFuncModel_nd : public FunctionModel_nd<t_real>
{
public:
	virtual ~FitterFuncModel_nd() = default;

	virtual std::size_t GetDim() const = 0;

	virtual bool SetParams(const std::vector<t_real>& vecParams) = 0;
	virtual t_real operator()(const std::vector<t_real>& dx) const = 0;

	virtual FitterFuncModel_nd* copy() const = 0;
	virtual std::string print(bool bFillInSyms=true) const = 0;

	virtual const char* GetModelName() const = 0;


	friend std::ostream& operator<<(std::ostream& ostr, const FitterFuncModel_nd<t_real>& fkt)
	{
		ostr << fkt.print();
		return ostr;
	}
};



// ----------------------------------------------------------------------------
// interface using supplied functions

/**
 * iNumArgs also includes the "x" parameter to the function, m_vecVals does not
 */
template<class t_real, std::size_t iNumArgs, typename t_func>
class FitterLamFuncModel : public FitterFuncModel<t_real>
{
protected:
	t_func m_func;
	std::vector<t_real> m_vecVals;
	bool m_bSeparateFreeParam = 1;	// separate "x" from parameters (for fitter)

public:
	FitterLamFuncModel(t_func func, bool bSeparateX=1)
		: m_func(func), m_bSeparateFreeParam(bSeparateX)
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
		return tl::call<iNumArgs, t_func, t_real, std::vector>(m_func, *pvecVals);
	}

	virtual FitterLamFuncModel* copy() const override
	{
		FitterLamFuncModel<t_real, iNumArgs, t_func>* pMod =
			new FitterLamFuncModel<t_real, iNumArgs, t_func>(m_func);

		pMod->m_vecVals = this->m_vecVals;
		pMod->m_bSeparateFreeParam = this->m_bSeparateFreeParam;

		return pMod;
	}

	virtual std::string print(bool bFillInSyms=true) const override { return ""; }
	virtual const char* GetModelName() const override { return "FitterLamFuncModel"; }

	virtual std::vector<std::string> GetParamNames() const override { return std::vector<std::string>(); }
	virtual std::vector<t_real> GetParamValues() const override { return m_vecVals; }
	virtual std::vector<t_real> GetParamErrors() const override { return std::vector<t_real>(); }
};


}

#endif
