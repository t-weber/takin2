/**
 * probability distributions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_DISTR_H__
#define __TLIBS_DISTR_H__

#include <vector>
#include <array>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/logistic.hpp>

namespace tl {
namespace m = boost::math;


enum class DistrType
{
	// continuous
	NORMAL,
	LOGNORMAL,
	CAUCHY,
	CHI2,
	STUDENT,
	FISHER,
	EXP,
	BETA,
	GAMMA,
	LOGISTIC,

	// discrete
	BINOMIAL,
	HYPERGEOMETRIC,
	POISSON,

	NONE,
};


// ----------------------------------------------------------------------------
/**
 * collect properties of the various distributions
 */
template<class t_real, class t_distr, class=void> struct distr_traits {};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::normal_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr bool bIsDiscrete = 0;
	static constexpr DistrType distr_type = DistrType::NORMAL;
	static constexpr const char* pcName = "Normal";
	static constexpr const char* pcParam1 = "mu";
	static constexpr const char* pcParam2 = "sigma";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::lognormal_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr bool bIsDiscrete = 0;
	static constexpr DistrType distr_type = DistrType::LOGNORMAL;
	static constexpr const char* pcName = "Log-Normal";
	static constexpr const char* pcParam1 = "mu";
	static constexpr const char* pcParam2 = "sigma";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::cauchy_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr bool bIsDiscrete = 0;
	static constexpr DistrType distr_type = DistrType::CAUCHY;
	static constexpr const char* pcName = "Cauchy";
	static constexpr const char* pcParam1 = "mu";
	static constexpr const char* pcParam2 = "sigma";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::poisson_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 1;
	static constexpr bool bIsDiscrete = 1;
	static constexpr DistrType distr_type = DistrType::POISSON;
	static constexpr const char* pcName = "Poisson";
	static constexpr const char* pcParam1 = "lambda";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::binomial_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr bool bIsDiscrete = 1;
	static constexpr DistrType distr_type = DistrType::BINOMIAL;
	static constexpr const char* pcName = "Binomial";
	static constexpr const char* pcParam1 = "n";
	static constexpr const char* pcParam2 = "p";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::hypergeometric_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 3;
	static constexpr bool bIsDiscrete = 1;
	static constexpr DistrType distr_type = DistrType::HYPERGEOMETRIC;
	static constexpr const char* pcName = "Hypergeometric";
	static constexpr const char* pcParam1 = "r";
	static constexpr const char* pcParam2 = "n";
	static constexpr const char* pcParam3 = "N";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::chi_squared_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 1;
	static constexpr bool bIsDiscrete = 0;
	static constexpr DistrType distr_type = DistrType::CHI2;
	static constexpr const char* pcName = "Chi^2";
	static constexpr const char* pcParam1 = "dof";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::students_t_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 1;
	static constexpr bool bIsDiscrete = 0;
	static constexpr DistrType distr_type = DistrType::STUDENT;
	static constexpr const char* pcName = "Student";
	static constexpr const char* pcParam1 = "dof";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::fisher_f_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr bool bIsDiscrete = 0;
	static constexpr DistrType distr_type = DistrType::FISHER;
	static constexpr const char* pcName = "Fisher";
	static constexpr const char* pcParam1 = "dof1";
	static constexpr const char* pcParam2 = "dof2";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::exponential_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 1;
	static constexpr bool bIsDiscrete = 0;
	static constexpr DistrType distr_type = DistrType::EXP;
	static constexpr const char* pcName = "Exponential";
	static constexpr const char* pcParam1 = "lambda";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::beta_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr bool bIsDiscrete = 0;
	static constexpr DistrType distr_type = DistrType::BETA;
	static constexpr const char* pcName = "Beta";
	static constexpr const char* pcParam1 = "alpha";
	static constexpr const char* pcParam2 = "beta";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::gamma_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr bool bIsDiscrete = 0;
	static constexpr DistrType distr_type = DistrType::GAMMA;
	static constexpr const char* pcName = "Gamma";
	static constexpr const char* pcParam1 = "mu";
	static constexpr const char* pcParam2 = "sigma";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::logistic_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr bool bIsDiscrete = 0;
	static constexpr DistrType distr_type = DistrType::LOGISTIC;
	static constexpr const char* pcName = "Logistic";
	static constexpr const char* pcParam1 = "mu";
	static constexpr const char* pcParam2 = "sigma";
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
/**
 * sort parameter names into runtime and static containers
 */
template<class t_distr_traits, class=void> struct _distr_params {};

template<class t_distr_traits>
struct _distr_params<t_distr_traits,
	typename std::enable_if<t_distr_traits::iNumArgs==1>::type>
{
	static constexpr std::array<const char*, t_distr_traits::iNumArgs> get_arr()
	{
		return std::array<const char*, t_distr_traits::iNumArgs>
		({ t_distr_traits::pcParam1 });
	}

	template<template<class...> class t_vec=std::vector>
	static t_vec<const char*> get_vec()
	{
		return t_vec<const char*>({ t_distr_traits::pcParam1 });
	}
};

template<class t_distr_traits>
struct _distr_params<t_distr_traits,
	typename std::enable_if<t_distr_traits::iNumArgs==2>::type>
{
	static constexpr std::array<const char*, t_distr_traits::iNumArgs> get_arr()
	{
		return std::array<const char*, t_distr_traits::iNumArgs>
		({ t_distr_traits::pcParam1, t_distr_traits::pcParam2 });
	}

	template<template<class...> class t_vec=std::vector>
	static t_vec<const char*> get_vec()
	{
		return t_vec<const char*>({ t_distr_traits::pcParam1, t_distr_traits::pcParam2 });
	}
};

template<class t_distr_traits>
struct _distr_params<t_distr_traits,
	typename std::enable_if<t_distr_traits::iNumArgs==3>::type>
{
	static constexpr std::array<const char*, t_distr_traits::iNumArgs> get_arr()
	{
		return std::array<const char*, t_distr_traits::iNumArgs>
		({ t_distr_traits::pcParam1, t_distr_traits::pcParam2, t_distr_traits::pcParam3 });
	}

	template<template<class...> class t_vec=std::vector>
	static t_vec<const char*> get_vec()
	{
		return t_vec<const char*>({ t_distr_traits::pcParam1, t_distr_traits::pcParam2,
			t_distr_traits::pcParam3 });
	}
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
/**
 * probability distribution base class
 */
template<class t_real=double>
class DistrBase
{
public:
	virtual t_real pdf(t_real x) const = 0;
	virtual t_real cdf(t_real x) const = 0;	// cdf(x) == P(X <= x)
	virtual t_real cdf_inv(t_real p) const = 0;

	virtual t_real operator()(t_real x) const { return pdf(x); }
};


/**
 * class for specific distributions
 */
template<class t_distr, class t_real=typename t_distr::value_type,
	std::size_t iParams=distr_traits<t_real, t_distr>::iNumArgs>
class Distr : public DistrBase<t_real>
{
public:
	using value_type = t_real;
	using distr_type = t_distr;
	using traits_type = distr_traits<t_real, t_distr>;

protected:
	t_distr distr;

public:
	template<std::size_t _iParams=iParams>
	Distr(t_real dParam1,
		typename std::enable_if<_iParams==1, void>::type* =nullptr)
		: distr(dParam1)
	{}

	template<std::size_t _iParams=iParams>
	Distr(t_real dParam1, t_real dParam2,
		typename std::enable_if<_iParams==2, void>::type* =nullptr)
		: distr(dParam1, dParam2)
	{}

	template<std::size_t _iParams=iParams>
	Distr(t_real dParam1, t_real dParam2, t_real dParam3,
		typename std::enable_if<_iParams==3, void>::type* =nullptr)
		: distr(dParam1, dParam2, dParam3)
	{}


	static constexpr const char* GetName()
	{ return traits_type::pcName; }

	static constexpr std::size_t GetNumParams()
	{ return traits_type::iNumArgs; }

	static constexpr std::array<const char*, GetNumParams()> GetParamNames()
	{ return _distr_params<traits_type>::get_arr(); }
	static std::vector<const char*> GetParamNamesVec()
	{ return _distr_params<traits_type>::get_vec(); }


	virtual t_real pdf(t_real x) const override
	{
		if(traits_type::bIsDiscrete) x = std::round(x);
		return m::pdf(distr, x);
	}

	virtual t_real cdf(t_real x) const override
	{
		if(traits_type::bIsDiscrete) x = std::round(x);
		return m::cdf(distr, x);
	}

	virtual t_real cdf_inv(t_real p) const override
	{
		return m::quantile(distr, p);
	}
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
/**
 * typedefs for specific distributions
 */
// continuous
template<class t_real> using t_normal_dist = Distr<m::normal_distribution<t_real>>;
template<class t_real> using t_lognormal_dist = Distr<m::lognormal_distribution<t_real>>;
template<class t_real> using t_cauchy_dist = Distr<m::cauchy_distribution<t_real>>;
template<class t_real> using t_chi2_dist = Distr<m::chi_squared_distribution<t_real>>;
template<class t_real> using t_student_dist = Distr<m::students_t_distribution<t_real>>;
template<class t_real> using t_fisher_dist = Distr<m::fisher_f_distribution<t_real>>;
template<class t_real> using t_exp_dist = Distr<m::exponential_distribution<t_real>>;
template<class t_real> using t_beta_dist = Distr<m::beta_distribution<t_real>>;
template<class t_real> using t_gamma_dist = Distr<m::gamma_distribution<t_real>>;
template<class t_real> using t_logistic_dist = Distr<m::logistic_distribution<t_real>>;

// discrete
template<class t_real> using t_poisson_dist = Distr<m::poisson_distribution<t_real>>;
template<class t_real> using t_binomial_dist = Distr<m::binomial_distribution<t_real>>;
template<class t_real> using t_hypergeo_dist = Distr<m::hypergeometric_distribution<t_real>>;
// ----------------------------------------------------------------------------


}
#endif
