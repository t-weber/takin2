/**
 * spins
 * @author: Tobias Weber <tobias.weber@tum.de>
 * @date: 2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_SPIN_H__
#define __TLIBS_SPIN_H__

#include "../math/linalg.h"
#include <boost/math/special_functions/factorials.hpp>


namespace tl {

/**
 * spin matrices
 * @see e.g. (Arfken 2013), p. 110
 */
template<template<class...> class t_mat=ublas::matrix,
	template<class...> class t_vec=ublas::vector,
	class t_real = double>
t_vec<t_mat<std::complex<t_real>>> get_spin_matrices()
{
	t_vec<t_mat<std::complex<t_real>>> vec(3);
	const std::complex<t_real> i(0,1);

	vec[0] = make_mat<t_mat<std::complex<t_real>>>({{0,1}, {1,0}});
	vec[1] = make_mat<t_mat<std::complex<t_real>>>({{0,-i}, {i,0}});
	vec[2] = make_mat<t_mat<std::complex<t_real>>>({{1,0}, {0,-1}});

	return vec;
}


template<template<class...> class t_mat=ublas::matrix,
	template<class...> class t_vec=ublas::vector,
	class t_real = double>
t_vec<t_mat<std::complex<t_real>>> get_ladder_ops()
{
	t_vec<t_mat<std::complex<t_real>>> vecS = get_spin_matrices();

	t_vec<t_mat<std::complex<t_real>>> vec(2);
	const std::complex<t_real> i(0,1);

	vec[0] = vecS[0] + i*vecS[1];	// up
	vec[1] = vecS[0] - i*vecS[1];	// down

	return vec;
}


template<class t_mat=ublas::matrix<double>>
t_mat commutator(const t_mat& A, const t_mat& B)
{
	t_mat AB = prod_mm(A, B);
	t_mat BA = prod_mm(B, A);
	return AB - BA;
}


/**
 * spin rotation in SU(2)
 * @see e.g. (Arfken 2013), p. 851
 */
template<template<class...> class t_mat = ublas::matrix, class t_real = double>
t_mat<std::complex<t_real>> rot_spin(int iComp, t_real dAngle)
{
	const auto vecS = get_spin_matrices<t_mat, ublas::vector, t_real>();
	const auto matI = unit_m<t_mat<std::complex<t_real>>>(2);
	const std::complex<t_real> I(0,1);

	t_mat<std::complex<t_real>> mat =
		std::complex<t_real>(std::cos(t_real(0.5)*dAngle)) * matI +
		std::complex<t_real>(std::sin(t_real(0.5)*dAngle)) * I*vecS[iComp];
	return mat;
}


/**
 * CG coefficients
 * @see (Arfken 2013), p. 790 for the formula
 * 
 * e.g. two e- spins: s1 = s2 = 0.5, ms[1,2] = 0.5 (up) or -0.5 (down), S = 0 (sing.) or 1 (trip.)
 */
template<class T = double>
T CG_coeff(T S, T s1, T s2, T ms1, T ms2)
{
	T (*fak)(T) = [](T t) -> T { return boost::math::factorial<T>(t); };

	T tCG = fak(S + s1 - s2)*fak(S - s1 + s2)*fak(-S + s1 + s2);
	tCG *= (T(2)*S + T(1));
	tCG *= fak(S + ms1 + ms2) * fak(S - (ms1 + ms2));
	tCG *= fak(s1 + ms1) * fak(s1 - ms1);
	tCG *= fak(s2 + ms2) * fak(s2 - ms2);
	tCG /= fak(S + s1 + s2 + T(1));
	tCG = std::sqrt(tCG);

	auto k_fkt = [&](T k) -> T
	{
		T t = std::pow(T(-1), k);
		t /= fak(k);
		t /= fak(-S + s1 + s2 - k)*fak(S - s1 - ms2 + k)*fak(S - s2 + ms1 + k);
		t /= fak(s2 + ms2 - k)*fak(s1 - ms1 - k);
		return t;
	};

	auto k_minmax = [&]() -> std::pair<T,T>
	{
		T kmax = s1 - ms1;
		kmax = std::min(kmax, s2 + ms2);
		kmax = std::min(kmax, -S + s1 + s2);

		T kmin = -(S - s1 - ms2);
		kmin = std::max(kmin, -(S - s2 + ms1));
		kmin = std::max(kmin, T(0));

		return std::make_pair(kmin, kmax);
	};

	T kmin, kmax;
	std::tie(kmin, kmax) = k_minmax();
	T kfact = T(0);
	for(T k=kmin; k<=kmax; k+=T(1))
		kfact += k_fkt(k);
	tCG *= kfact;

	return tCG;
}

}

#endif
