/**
 * statistical functions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_STAT_H__
#define __TLIBS_STAT_H__

#include "linalg.h"
#include "distr.h"
#include "numint.h"

#include <boost/math/special_functions/binomial.hpp>
//#include <boost/numeric/ublas/io.hpp>


namespace tl {

/**
 * mean value <x>
 */
template<class vec_type>
typename vec_type::value_type mean_value(const vec_type& vec)
{
	using T = typename vec_type::value_type;
	if(vec.size()==0) return T(0);
	else if(vec.size()==1) return vec[0];

	T tMean = vec[0];
	for(std::size_t i=1; i<vec.size(); ++i)
		tMean += vec[i];
	tMean /= vec.size();

	return tMean;
}


/**
 * mean value of <x^2>
 */
template<class vec_type>
typename vec_type::value_type mean_square_value(const vec_type& vec)
{
	using T = typename vec_type::value_type;
	if(vec.size()==0) return T(0);
	else if(vec.size()==1) return vec[0];

	T tMean = vec[0]*vec[0];
	for(std::size_t i=1; i<vec.size(); ++i)
		tMean += vec[i]*vec[i];
	tMean /= vec.size();

	return tMean;
}


/**
 * mean value with given probability
 */
template<class vec_type_prob, class vec_type>
typename vec_type::value_type mean_value(const vec_type_prob& vecP, const vec_type& vec)
{
	typedef typename vec_type::value_type T;
	typedef typename vec_type_prob::value_type Tprob;
	std::size_t iSize = std::min(vecP.size(), vec.size());

	if(iSize==0) return T(0);

	T tMean = vecP[0]*vec[0];
	Tprob tProbTotal = vecP[0];
	for(std::size_t i=1; i<iSize; ++i)
	{
		tMean += vecP[i]*vec[i];
		tProbTotal += vecP[i];
	}
	tMean /= tProbTotal;

	return tMean;
}


/**
 * standard deviation of mean value, with correction factor
 * @see e.g.: https://en.wikipedia.org/wiki/Bessel%27s_correction
 */
template<class vec_type>
typename vec_type::value_type std_dev(const vec_type& vec, bool bCorr=1)
{
	typedef typename vec_type::value_type T;
	if(vec.size()<=1) return T(0);

	T tProb = T(vec.size());
	if(bCorr) tProb -= T(1);

	T tMean = mean_value(vec);
	T t = T(0);
	for(const T& tval : vec)
		t += (tval-tMean) * (tval-tMean);
	t /= tProb;

	return std::sqrt(t);
}


/**
 * standard deviation with given probability
 */
template<class vec_type_prob, class vec_type>
typename vec_type::value_type std_dev(const vec_type_prob& vecP, const vec_type& vec)
{
	typedef typename vec_type::value_type T;
	std::size_t iSize = std::min(vecP.size(), vec.size());
	if(iSize<=1) return T(0);

	T tMean = mean_value<vec_type_prob, vec_type>(vecP, vec);
	T t = T(0);
	T tProbTotal = T(0);

	for(std::size_t iIdx = 0; iIdx<iSize; ++iIdx)
	{
		t += (vec[iIdx]-tMean)*(vec[iIdx]-tMean) * vecP[iIdx];
		tProbTotal += vecP[iIdx];
	}
	t /= tProbTotal;

	return std::sqrt(t);
}


/**
 * entropy of a discrete distribution
 * S = - < log p_i >
 * @see e.g.: https://en.wikipedia.org/wiki/Entropy_(information_theory)
 */
template<class t_real=double, class t_func>
t_real entropy(const t_func& funcPdf, std::size_t iMax)
{
	t_real dS = 0;
	for(std::size_t iArg=0; iArg<iMax; ++iArg)
	{
		t_real dVal = funcPdf(t_real(iArg));
		if(!float_equal(dVal, t_real(0)))
			dS += dVal * std::log(dVal);
	}
	return -dS;
}


/**
 * entropy of a continuous distribution
 * S = - < log p(x_i) >
 * @see e.g.: https://en.wikipedia.org/wiki/Entropy_(information_theory)
 */
template<class t_real=double, class t_func>
t_real entropy(const t_func& funcPdf, t_real dXMin, t_real dXMax, std::size_t iSteps=128)
{
	std::function<t_real(t_real)> fktInt = [&funcPdf](t_real dX) -> t_real
	{
		t_real dVal = funcPdf(dX);
		if(!float_equal(dVal, t_real(0)))
			return dVal * std::log(dVal);
		return t_real(0);
	};

	t_real dS = numint_simpN(fktInt, dXMin, dXMax, iSteps);
	return -dS;
}


// -----------------------------------------------------------------------------


/**
 * Stirling's formula for log(n!)
 * @see e.g. https://en.wikipedia.org/wiki/Stirling%27s_approximation
 */
template<class t_real = double>
t_real log_nfac(t_real n)
{
	const t_real twopi = t_real(2) * get_pi<t_real>();
	return n*std::log(n) - n + std::log(twopi*t_real(n)) / t_real(2);
}


/**
 * combinatorics
 * @see e.g.: https://de.wikipedia.org/wiki/Abz%C3%A4hlende_Kombinatorik
 */
template<class t_real = double, class t_uint = unsigned>
t_real combinatorics(t_uint n, t_uint k, bool bOrdered, bool bRepetition)
{
	t_real tVal = t_real(0);

	if(bOrdered)	// variation
	{
		if(bRepetition)	// repetition of particles
		{
			tVal = std::pow(t_real(n), t_real(k));
		}
		else	// no repetition of particles
		{
			t_real binom = boost::math::binomial_coefficient<t_real>(n, k);
			tVal = binom * boost::math::factorial<t_real>(k);
		}
	}
	else			// combination
	{
		if(bRepetition)
			tVal = boost::math::binomial_coefficient<t_real>(n+k-1, k);
		else
			tVal = boost::math::binomial_coefficient<t_real>(n, k);
	}

	return tVal;
}


/**
 * possibilities to distribute particles onto niveaus
 * @see e.g.: https://de.wikipedia.org/wiki/Abz%C3%A4hlende_Kombinatorik
 */
template<class t_real = double, class t_uint = unsigned>
t_real particles_in_niveaus(t_uint iPart, t_uint iNiv, bool bDistinct, bool bOnePerNiveau)
{
	t_real tCnt = t_real(0);

	if(bDistinct)	// classical
	{
		if(bOnePerNiveau)
			tCnt = combinatorics<t_real, t_uint>(iNiv, iPart, true, false);
		else	// Boltzons
			tCnt = combinatorics<t_real, t_uint>(iNiv, iPart, true, true);
	}
	else	// qm
	{
		if(bOnePerNiveau)	// Fermions
			tCnt = combinatorics<t_real, t_uint>(iNiv, iPart, false, false);
		else	// Bosons
			tCnt = combinatorics<t_real, t_uint>(iNiv, iPart, false, true);
	}

	return tCnt;
}


template<class t_real = double, class t_uint = unsigned>
t_real bosons_in_niveaus(t_uint iPart, t_uint iNiv)
{ return particles_in_niveaus<t_real, t_uint>(iPart, iNiv, 0, 0); }


template<class t_real = double, class t_uint = unsigned>
t_real fermions_in_niveaus(t_uint iPart, t_uint iNiv)
{ return particles_in_niveaus<t_real, t_uint>(iPart, iNiv, 0, 1); }


template<class t_real = double, class t_uint = unsigned>
t_real boltzons_in_niveaus(t_uint iPart, t_uint iNiv)
{ return particles_in_niveaus<t_real, t_uint>(iPart, iNiv, 1, 0); }


// -----------------------------------------------------------------------------


/**
 * calculates the covariance and the correlation matrices
 * covariance: C_ij = cov(X_i, X_j) = < (X_i - <X_i>) * (X_j - <X_j>) >
 * correlation: K_ij = C_ij / (sigma_i sigma_j)
 * @see e.g.: http://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm
 * @see e.g.: (Arfken 2013) p. 1142
 */
template<typename T=double>
std::tuple<ublas::matrix<T>, ublas::matrix<T>>
covariance(const std::vector<ublas::vector<T>>& vecVals, const std::vector<T>* pProb = 0)
{
	using t_mat = ublas::matrix<T>;
	using t_vecvec = typename std::remove_reference<decltype(vecVals)>::type;
	using t_innervec_org = decltype(vecVals[0]);
	using t_innervec = typename std::remove_const<
		typename std::remove_reference<t_innervec_org>::type>::type;

	if(vecVals.size() == 0) return std::make_tuple(t_mat(), t_mat());

	// mean vector <X_i>
	t_innervec vecMean;
	if(pProb)
		vecMean = mean_value<std::vector<T>, t_vecvec>(*pProb, vecVals);
	else
		vecMean = mean_value<t_vecvec>(vecVals);

	t_mat matCov = zero_m<t_mat>(vecVals[0].size(), vecVals[0].size());
	T tSum = T(0);
	const std::size_t N = vecVals.size();

	for(std::size_t i=0; i<N; ++i)
	{
		T tprob = T(1);

		// X_i - <X_i>
		t_innervec vec = vecVals[i] - vecMean;

		// matrix elements, AA^t
		t_mat matOuter = outer(vec, vec);

		// probabilities for final averaging, <...>
		if(pProb)
		{
			tprob = (*pProb)[i];
			matOuter *= tprob;
		}

		matCov += matOuter;
		tSum += tprob;
	}

	// average, sometimes defined as C /= (N-1)
	matCov /= tSum /*-T(1)*/;


	// --------------------------------------------------------------------------------
	// correlation matrix
	t_innervec vecVar = diag_vec(matCov);
	t_innervec vecStdDev(vecVar.size());

	std::transform(vecVar.begin(), vecVar.end(), vecStdDev.begin(),
		[](typename t_innervec::value_type d) -> typename t_innervec::value_type
		{ return std::sqrt(d); });

	t_mat matStdDev = outer(vecStdDev, vecStdDev);
	t_mat matCorr = ublas::element_div(matCov, matStdDev);
	// --------------------------------------------------------------------------------

	return std::make_tuple(matCov, matCorr);
}


// -----------------------------------------------------------------------------


/**
 * calculates chi^2 distance of a function model to data points
 * chi^2 = sum( (y_i - f(x_i))^2 / sigma_i^2 )
 * @see e.g.: (Arfken 2013), p. 1170
 */
template<class T, class t_func, class t_iter_dat=T*>
T chi2(const t_func& func, std::size_t N,
	const t_iter_dat x, const t_iter_dat y, const t_iter_dat dy)
{
	using t_dat = typename std::remove_pointer<t_iter_dat>::type;
	T tchi2 = T(0);

	for(std::size_t i=0; i<N; ++i)
	{
		T td = T(y[i]) - func(T(x[i]));
		T tdy = dy ? T(dy[i]) : T(0.1*td);	// 10% error if none given

		if(std::abs(tdy) < std::numeric_limits<t_dat>::min())
			tdy = std::numeric_limits<t_dat>::min();

		T tchi = T(td) / T(tdy);
		tchi2 += tchi*tchi;
	}

	return tchi2;
}

template<class t_vec, class t_func>
typename t_vec::value_type chi2(const t_func& func,
	const t_vec& x, const t_vec& y, const t_vec& dy)
{
	using T = typename t_vec::value_type;
	return chi2<T, t_func, T*>(func, x.size(), x.data(), y.data(),
		dy.size() ? dy.data() : nullptr);
}


/**
 * chi^2 which doesn't use an x value, but an index instead: y[idx] - func(idx)
 */
template<class T, class t_func, class t_iter_dat=T*>
T chi2_idx(const t_func& func, std::size_t N, const t_iter_dat y, const t_iter_dat dy)
{
	using t_dat = typename std::remove_pointer<t_iter_dat>::type;
	T tchi2 = T(0);

	for(std::size_t i=0; i<N; ++i)
	{
		T td = T(y[i]) - func(i);
		T tdy = dy ? T(dy[i]) : T(0.1*td);	// 10% error if none given

		if(std::abs(tdy) < std::numeric_limits<t_dat>::min())
			tdy = std::numeric_limits<t_dat>::min();

		T tchi = T(td) / T(tdy);
		tchi2 += tchi*tchi;
	}

	return tchi2;
}


/**
 * direct chi^2 calculation with a model array instead of a model function
 */
template<class T, class t_iter_dat=T*>
T chi2_direct(std::size_t N, const t_iter_dat func_y, const t_iter_dat y, const t_iter_dat dy)
{
	using t_dat = typename std::remove_pointer<t_iter_dat>::type;
	T tchi2 = T(0);

	for(std::size_t i=0; i<N; ++i)
	{
		T td = T(y[i]) - T(func_y[i]);
		T tdy = dy ? T(dy[i]) : T(0.1*td);	// 10% error if none given

		if(std::abs(tdy) < std::numeric_limits<t_dat>::min())
			tdy = std::numeric_limits<t_dat>::min();

		T tchi = T(td) / T(tdy);
		tchi2 += tchi*tchi;
	}

	return tchi2;
}


/**
 * multi-dimensional chi^2 function
 */
template<class T, class T_dat, class t_func, template<class...> class t_vec=std::vector>
T chi2_nd(const t_func& func,
	const t_vec<t_vec<T_dat>>& vecvecX, const t_vec<T_dat>& vecY, const t_vec<T_dat>& vecDY)
{
	T tchi2 = T(0);

	for(std::size_t i=0; i<vecvecX.size(); ++i)
	{
		T td = T(vecY[i]) - func(vecvecX[i]);
		T tdy = vecDY[i];

		if(std::abs(tdy) < std::numeric_limits<T_dat>::min())
			tdy = std::numeric_limits<T_dat>::min();

		T tchi = T(td) / T(tdy);
		tchi2 += tchi*tchi;
	}

	return tchi2;
}


// -----------------------------------------------------------------------------


/**
 * Confidence interval of array data mean using t-distribution
 * @see e.g.: (Arfken 2013), pp. 1176ff
 */
template<class t_real = double, class t_vec = std::vector<t_real>>
std::tuple<t_real, t_real, t_real>	// [mean, stddev, confidence]
confidence(const t_vec& vec, t_real dProb)
{
	t_real dMean = mean_value(vec);
	t_real dStd = std_dev(vec, true);

	t_real dDof = t_real(vec.size()-1);
	tl::t_student_dist<t_real> t(dDof);
	t_real dConf = t.cdf_inv(0.5 + dProb/2.);
	dConf *= dStd / std::sqrt(dDof);

	return std::make_tuple(dMean, dStd, dConf);
}

}

#endif
