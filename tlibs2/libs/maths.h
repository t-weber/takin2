/**
 * tlibs2
 * (container-agnostic) math library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * @note The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 * @note Additional functions forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @note Further functions and updates forked on 1-Feb-2021 and 19-Apr-2021 from my privately developed "geo" and "misc" projects (https://github.com/t-weber/geo and https://github.com/t-weber/misc).
 * @note Additional functions forked on 6-Feb-2022 from my privately developed "mathlibs" project (https://github.com/t-weber/mathlibs).
 *
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "magtools", "geo", "misc", and "mathlibs" projects
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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

#ifndef __TLIBS2_CXX20_MATH_ALGOS_H__
#define __TLIBS2_CXX20_MATH_ALGOS_H__

//#define USE_LAPACK 1
#define __TLIBS2_QR_METHOD 0

#include <cstddef>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <complex>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <initializer_list>
#include <limits>
#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <random>
#include <utility>
#include <optional>
#include <memory>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#if __has_include(<numbers>)
	#include <numbers>

	#define __TLIBS2_USE_NUMBERS__
#endif

#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/algorithm/minmax_element.hpp>

#include "log.h"
#include "str.h"
#include "bits.h"
#include "traits.h"

#if __has_include(<lapacke.h>) && USE_LAPACK
	extern "C"
	{
		#define lapack_complex_double std::complex<double>
		#define lapack_complex_double_real(z) (z.real())
		#define lapack_complex_double_imag(z) (z.imag())

		#define lapack_complex_float std::complex<float>
		#define lapack_complex_float_real(z) (z.real())
		#define lapack_complex_float_imag(z) (z.imag())

		#include <lapacke.h>
	}

	#define __TLIBS2_USE_LAPACK__
#else
	#pragma message("tlibs2: Disabling Lapack(e) library (not found).")
#endif

#if __has_include(<Faddeeva.hh>)
	#include <Faddeeva.hh>
	using t_real_fadd = double;

	#define __TLIBS2_USE_FADDEEVA__
#else
	#pragma message("tlibs2: Disabling Faddeeva library (not found).")
#endif

#if __has_include(<Qhull.h>)
	#include <Qhull.h>
	#include <QhullFacetList.h>
	#include <QhullVertexSet.h>

	#define __TLIBS2_USE_QHULL__
#else
	#pragma message("tlibs2: Disabling QHull library (not found).")
#endif

// separator tokens
#define TL2_COLSEP ';'
#define TL2_ROWSEP '|'


namespace tl2 {
// ----------------------------------------------------------------------------
// forward declarations
// ----------------------------------------------------------------------------
template<class t_mat>
t_mat prod(const t_mat& mat1, const t_mat& mat2, bool assert_sizes=true)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>;

template<class t_mat>
t_mat zero(std::size_t N1, std::size_t N2)
requires is_basic_mat<t_mat>;

template<class t_vec>
t_vec zero(std::size_t N=0)
requires is_basic_vec<t_vec>;

template<class t_mat, class t_vec>
std::tuple<bool, t_mat, t_mat> qr(const t_mat& mat)
requires is_mat<t_mat> && is_vec<t_vec>;

template<class t_mat>
std::tuple<t_mat, bool> inv(const t_mat& mat)
requires is_mat<t_mat>;

template<class t_mat>
typename t_mat::value_type det(const t_mat& mat)
requires is_mat<t_mat>;
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// helpers and constants
// ----------------------------------------------------------------------------

// constants
#ifdef __TLIBS2_USE_NUMBERS__
	template<typename T=double> constexpr T golden{std::numbers::phi_v<T>};
	template<typename T=double> constexpr T pi{std::numbers::pi_v<T>};
#else
	// see: https://en.wikipedia.org/wiki/Golden_ratio
	template<typename T=double> constexpr T golden{1.618033988749895};
	template<typename T=double> constexpr T pi{M_PI};
#endif

template<typename INT=int> bool is_even(INT i) { return (i%2 == 0); }
template<typename INT=int> bool is_odd(INT i) { return !is_even<INT>(i); }

template<class T=double> constexpr T r2d(T rad) { return rad/pi<T>*T(180); }	// rad -> deg
template<class T=double> constexpr T d2r(T deg) { return deg/T(180)*pi<T>; }	// deg -> rad
template<class T=double> constexpr T r2m(T rad) { return rad/pi<T>*T(180*60); }	// rad -> min
template<class T=double> constexpr T m2r(T min) { return min/T(180*60)*pi<T>; }	// min -> rad

/**
 * Gaussian around 0: f(x) = exp(-1/2 * (x/sig)^2)
 * at hwhm: f(x_hwhm) = 1/2
 *          exp(-1/2 * (x_hwhm/sig)^2) = 1/2
 *          -1/2 * (x_hwhm/sig)^2 = ln(1/2)
 *          (x_hwhm/sig)^2 = -2*ln(1/2)
 *          x_hwhm^2 = sig^2 * 2*ln(2)
 */
template<class T=double> static constexpr T SIGMA2FWHM = T(2)*std::sqrt(T(2)*std::log(T(2)));
template<class T=double> static constexpr T SIGMA2HWHM = std::sqrt(T(2)*std::log(T(2)));
template<class T=double> static constexpr T FWHM2SIGMA = T(1)/SIGMA2FWHM<T>;
template<class T=double> static constexpr T HWHM2SIGMA = T(1)/SIGMA2HWHM<T>;


template<typename T> T sign(T t)
{
	if(t<0.) return -T(1);
	return T(1);
}


template<typename T> T cot(T t)
{
	return std::tan(T(0.5)*pi<T> - t);
}


template<typename T> T coth(T t)
{
	return T(1) / std::tanh(t);
}


template<typename T=double>
T log(T tbase, T tval)
{
	return T(std::log(tval)/std::log(tbase));
}


template<typename T=double>
T nextpow(T tbase, T tval)
{
	return T(std::pow(tbase, std::ceil(log(tbase, tval))));
}


/**
 * unsigned angle between two vectors
 * <q1|q2> / (|q1| |q2|) = cos(alpha)
 */
template<class t_vec>
typename t_vec::value_type angle_unsigned(const t_vec& q1, const t_vec& q2)
requires is_basic_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	if(q1.size() != q2.size())
		return t_real(0);

	t_real dot = t_real(0);
	t_real len1 = t_real(0);
	t_real len2 = t_real(0);

	for(std::size_t i=0; i<q1.size(); ++i)
	{
		dot += q1[i]*q2[i];

		len1 += q1[i]*q1[i];
		len2 += q2[i]*q2[i];
	}

	len1 = std::sqrt(len1);
	len2 = std::sqrt(len2);

	dot /= len1;
	dot /= len2;

	return std::acos(dot);
}


/**
 * unsigned angle between two quaternions
 * <q1|q2> / (|q1| |q2|) = cos(alpha)
 */
template<class t_quat>
typename t_quat::value_type angle_unsigned(const t_quat& q1, const t_quat& q2)
requires is_quat<t_quat>
{
	using t_real = typename t_quat::value_type;

	t_real dot = q1.R_component_1() * q2.R_component_1() +
		q1.R_component_2() * q2.R_component_2() +
		q1.R_component_3() * q2.R_component_3() +
		q1.R_component_4() * q2.R_component_4();

	t_real len1 = q1.R_component_1() * q1.R_component_1() +
		q1.R_component_2() * q1.R_component_2() +
		q1.R_component_3() * q1.R_component_3() +
		q1.R_component_4() * q1.R_component_4();

	t_real len2 = q2.R_component_1() * q2.R_component_1() +
		q2.R_component_2() * q2.R_component_2() +
		q2.R_component_3() * q2.R_component_3() +
		q2.R_component_4() * q2.R_component_4();

	len1 = std::sqrt(len1);
	len2 = std::sqrt(len2);

	dot /= len1;
	dot /= len2;

	return std::acos(dot);
}


/**
 * slerp
 * @see K. Shoemake, "Animating rotation with quaternion curves", http://dx.doi.org/10.1145/325334.325242
 * @see (Desktop Bronstein 2008), formula 4.207
 * @see (Bronstein 2008), p. 306, formula 4.155
 */
template<class T>
T slerp(const T& q1, const T& q2, typename T::value_type t)
{
	using t_real = typename T::value_type;
	t_real angle = angle_unsigned<T>(q1, q2);

	T q = std::sin((t_real(1)-t)*angle)/std::sin(angle) * q1 +
	std::sin(t*angle)/std::sin(angle) * q2;

	return q;
}


/**
 * x = 0..1, y = 0..1
 * @see https://en.wikipedia.org/wiki/Bilinear_interpolation
 */
template<typename T=double>
T bilinear_interp(T x0y0, T x1y0, T x0y1, T x1y1, T x, T y)
{
	T top = std::lerp(x0y1, x1y1, x);
	T bottom = std::lerp(x0y0, x1y0, x);

	return std::lerp(bottom, top, y);
}


template<typename T=double, typename REAL=double,
template<class...> class t_vec = std::vector>
t_vec<T> linspace(const T& tmin, const T& tmax, std::size_t iNum)
{
	t_vec<T> vec;
	vec.reserve(iNum);

	if(iNum == 1)
	{
		// if just one point is requested, use the lower limit
		vec.push_back(tmin);
		return vec;
	}

	for(std::size_t i=0; i<iNum; ++i)
		vec.push_back(std::lerp(tmin, tmax, REAL(i)/REAL(iNum-1)));
	return vec;
}


template<typename T=double, typename REAL=double,
template<class...> class t_vec = std::vector>
t_vec<T> logspace(const T& tmin, const T& tmax, std::size_t iNum, T tBase=T(10))
{
	t_vec<T> vec = linspace<T, REAL>(tmin, tmax, iNum);

	for(T& t : vec)
		t = std::pow(tBase, t);
	return vec;
}


template<typename T>
T clamp(T t, T min, T max)
{
	if(t < min) t = min;
	if(t > max) t = max;

	return t;
}


template<class T>
bool is_in_range(T val, T centre, T pm)
{
	pm = std::abs(pm);

	if(val < centre-pm) return false;
	if(val > centre+pm) return false;
	return true;
}


/**
 * point contained in linear range?
 */
template<class T = double>
bool is_in_linear_range(T dStart, T dStop, T dPoint)
{
	if(dStop < dStart)
		std::swap(dStart, dStop);

	return (dPoint >= dStart) && (dPoint <= dStop);
}


/**
 * angle contained in angular range?
 */
template<class T = double>
bool is_in_angular_range(T dStart, T dRange, T dAngle)
{
	if(dStart < T(0)) dStart += T(2)*pi<T>;
	if(dAngle < T(0)) dAngle += T(2)*pi<T>;

	dStart = std::fmod(dStart, T(2)*pi<T>);
	dAngle = std::fmod(dAngle, T(2)*pi<T>);

	T dStop = dStart + dRange;


	// if the end point is contained in the circular range
	if(dStop < T(2)*pi<T>)
	{
		return is_in_linear_range<T>(dStart, dStop, dAngle);
	}
	// else end point wraps around
	else
	{
		return is_in_linear_range<T>(dStart, T(2)*pi<T>, dAngle) ||
			is_in_linear_range<T>(T(0), dRange-(T(2)*pi<T>-dStart), dAngle);
	}
}


/**
 * converts a string to a scalar value
 */
template<class t_scalar=double, class t_str=std::string>
t_scalar stoval(const t_str& str)
{
	if constexpr(std::is_same_v<t_scalar, float>)
		return std::stof(str);
	else if constexpr(std::is_same_v<t_scalar, double>)
		return std::stod(str);
	else if constexpr(std::is_same_v<t_scalar, long double>)
		return std::stold(str);
	else if constexpr(std::is_same_v<t_scalar, int>)
		return std::stoi(str);
//	else if constexpr(std::is_same_v<t_scalar, unsigned int>)
//		return std::stoui(str);
	else if constexpr(std::is_same_v<t_scalar, long>)
		return std::stol(str);
	else if constexpr(std::is_same_v<t_scalar, unsigned long>)
		return std::stoul(str);
	else if constexpr(std::is_same_v<t_scalar, long long>)
		return std::stoll(str);
	else if constexpr(std::is_same_v<t_scalar, unsigned long long>)
		return std::stoull(str);
	else
	{
		t_scalar val{};
		std::istringstream{str} >> val;
		return val;
	}
}


/**
 * get a random number in the given range
 */
template<class t_num>
t_num get_rand(t_num min=1, t_num max=-1)
{
	static std::mt19937 rng{std::random_device{}()};

	if(max <= min)
	{
		min = std::numeric_limits<t_num>::lowest() / 10.;
		max = std::numeric_limits<t_num>::max() / 10.;
	}

	if constexpr(std::is_integral_v<t_num>)
		return std::uniform_int_distribution<t_num>(min, max)(rng);
	else
		return std::uniform_real_distribution<t_num>(min, max)(rng);
}
// ----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// peak model functions
// -----------------------------------------------------------------------------
/**
 * gaussian
 * @see https://en.wikipedia.org/wiki/Gaussian_function
 */
template<class T = double>
T gauss_model(T x, T x0, T sigma, T amp, T offs)
{
	T norm = T(1)/(std::sqrt(T(2)*pi<T>) * sigma);
	return amp * norm * std::exp(-0.5 * ((x-x0)/sigma)*((x-x0)/sigma)) + offs;
}


template<class T = double>
T gauss_model_amp(T x, T x0, T sigma, T amp, T offs)
{
	return amp * std::exp(-0.5 * ((x-x0)/sigma)*((x-x0)/sigma)) + offs;
}


template<class T = double>
T gauss_model_amp_slope(T x, T x0, T sigma, T amp, T offs, T slope)
{
	return amp * std::exp(-0.5 * ((x-x0)/sigma)*((x-x0)/sigma)) + (x-x0)*slope + offs;
}


/**
 * lorentzian
 * @see https://en.wikipedia.org/wiki/Cauchy_distribution
 */
template<class T = double>
T lorentz_model_amp(T x, T x0, T hwhm, T amp, T offs)
{
	return amp*hwhm*hwhm / ((x-x0)*(x-x0) + hwhm*hwhm) + offs;
}


template<class T = double>
T lorentz_model_amp_slope(T x, T x0, T hwhm, T amp, T offs, T slope)
{
	return amp*hwhm*hwhm / ((x-x0)*(x-x0) + hwhm*hwhm) + (x-x0)*slope + offs;
}


template<class T = double>
T parabola_model(T x, T x0, T amp, T offs)
{
	return amp*(x-x0)*(x-x0) + offs;
}


template<class T = double>
T parabola_model_slope(T x, T x0, T amp, T offs, T slope)
{
	return amp*(x-x0)*(x-x0) + (x-x0)*slope + offs;
}


// -----------------------------------------------------------------------------
template<class t_real_to, class t_real_from,
bool bIsEqu = std::is_same<t_real_from, t_real_to>::value>
struct complex_cast
{
	const std::complex<t_real_to>& operator()(const std::complex<t_real_from>& c) const
	{ return c; }
};

template<class t_real_to, class t_real_from>
struct complex_cast<t_real_to, t_real_from, 0>
{
	std::complex<t_real_to> operator()(const std::complex<t_real_from>& c) const
	{ return std::complex<t_real_to>(t_real_to(c.real()), t_real_to(c.imag())); }
};
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// Faddeeva function
// -----------------------------------------------------------------------------
#ifdef __TLIBS2_USE_FADDEEVA__

/**
 * Complex error function
 */
template<class T=double>
std::complex<T> erf(const std::complex<T>& z)
{
	complex_cast<t_real_fadd, T> cst;
	complex_cast<T, t_real_fadd> inv_cst;
	return inv_cst(::Faddeeva::erf(cst(z)));
}

/**
 * Complex complementary error function
 */
template<class T=double>
std::complex<T> erfc(const std::complex<T>& z)
{
	complex_cast<t_real_fadd, T> cst;
	complex_cast<T, t_real_fadd> inv_cst;
	return inv_cst(::Faddeeva::erfc(cst(z)));
}

/**
 * Faddeeva function
 * @see https://en.wikipedia.org/wiki/Faddeeva_function
 */
template<class T=double>
std::complex<T> faddeeva(const std::complex<T>& z)
{
	std::complex<T> i(0, 1.);
	return std::exp(-z*z) * erfc(-i*z);
}

/**
 * Voigt profile
 * @see https://en.wikipedia.org/wiki/Voigt_profile
 */
template<class T=double>
T voigt_model(T x, T x0, T sigma, T gamma, T amp, T offs)
{
	T norm = T(1)/(std::sqrt(T(2)*pi<T>) * sigma);
	std::complex<T> z = std::complex<T>(x-x0, gamma) / (sigma * std::sqrt(T(2)));

	return amp*norm * faddeeva<T>(z).real() + offs;
}

template<class T=double>
T voigt_model_amp(T x, T x0, T sigma, T gamma, T amp, T offs)
{
	std::complex<T> z = std::complex<T>(x-x0, gamma) / (sigma * std::sqrt(T(2)));
	return amp * faddeeva<T>(z).real() + offs;
}

template<class T=double>
T voigt_model_amp_slope(T x, T x0, T sigma, T gamma, T amp, T offs, T slope)
{
	std::complex<T> z = std::complex<T>(x-x0, gamma) / (sigma * std::sqrt(T(2)));
	return amp * faddeeva<T>(z).real() + (x-x0)*slope + offs;
}

#endif
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// physics-related functions
// -----------------------------------------------------------------------------

/**
 * wrapper for boost's Y function
 */
template<class T=double>
std::complex<T> Ylm(int l /*0..i*/, int m /*-l..l*/, T th /*0..pi*/, T ph /*0..2pi*/)
{
	return boost::math::spherical_harmonic<T,T>(l,m, th, ph);
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
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// coordinate trafos
// -----------------------------------------------------------------------------

/**
 * cartesian -> spherical
 * @see https://en.wikipedia.org/wiki/Spherical_coordinate_system
 */
template<class T = double>
std::tuple<T,T,T> cart_to_sph(T x, T y, T z)
{
	T rho = std::sqrt(x*x + y*y + z*z);
	T phi = std::atan2(y, x);
	T theta = std::acos(z/rho);

	return std::make_tuple(rho, phi, theta);
}


/**
 * spherical -> cartesian
 * @see https://en.wikipedia.org/wiki/Spherical_coordinate_system
 */
template<class T = double>
std::tuple<T,T,T> sph_to_cart(T rho, T phi, T theta)
{
	T x = rho * std::cos(phi)*std::sin(theta);
	T y = rho * std::sin(phi)*std::sin(theta);
	T z = rho * std::cos(theta);

	return std::make_tuple(x, y, z);
}


/**
 * cylindrical -> spherical
 * @see https://en.wikipedia.org/wiki/Spherical_coordinate_system
 */
template<class T = double>
std::tuple<T,T,T> cyl_to_sph(T rho_cyl, T phi_cyl, T z_cyl)
{
	T rho = std::sqrt(rho_cyl*rho_cyl + z_cyl*z_cyl);
	T theta = std::acos(z_cyl/rho);

	return std::make_tuple(rho, phi_cyl, theta);
}


/**
 * spherical -> cylindrical
 * @see https://en.wikipedia.org/wiki/Spherical_coordinate_system
 */
template<class T = double>
std::tuple<T,T,T> sph_to_cyl(T rho_sph, T phi_sph, T theta_sph)
{
	T rho = rho_sph * std::sin(theta_sph);
	T z = rho_sph * std::cos(theta_sph);

	return std::make_tuple(rho, phi_sph, z);
}


/**
 * cylindrical -> cartesian
 * @see https://en.wikipedia.org/wiki/Cylindrical_coordinate_system
 */
template<class T = double>
std::tuple<T,T,T> cyl_to_cart(T rho, T phi, T z)
{
	T x = rho * std::cos(phi);
	T y = rho * std::sin(phi);

	return std::make_tuple(x, y, z);
}


/**
 * cartesian -> cylindrical
 * @see https://en.wikipedia.org/wiki/Cylindrical_coordinate_system
 */
template<class T = double>
std::tuple<T,T,T> cart_to_cyl(T x, T y, T z)
{
	T rho = std::sqrt(x*x + y*y);
	T phi = std::atan2(y, x);

	return std::make_tuple(rho, phi, z);
}



template<class T = double>
std::tuple<T,T> crys_to_sph(T twophi_crys, T twotheta_crys)
{
	// converts the out-of-plane scattering angle '2theta' to the spherical theta
	T theta_sph = pi<T>/T(2) - twotheta_crys;
	// converts in-plane scattering angle '2phi' to the spherical phi
	T phi_sph = twophi_crys - pi<T>/T(2);

	return std::make_tuple(phi_sph, theta_sph);
}


template<class T = double>
std::tuple<T,T> sph_to_crys(T phi, T theta)
{
	return crys_to_sph<T>(phi, theta);
}


/**
 * gnomonic projection (similar to perspective projection with fov=90°)
 * @return [x,y]
 * @see see http://mathworld.wolfram.com/GnomonicProjection.html
 */
template<class T = double>
std::tuple<T,T> gnomonic_proj(T twophi_crys, T twotheta_crys)
{
	T x = -std::tan(twophi_crys);
	T y = std::tan(twotheta_crys) / std::cos(twophi_crys);

	return std::make_tuple(x, y);
}

/**
 * stereographic projection
 * @return [x,y]
 * @see http://mathworld.wolfram.com/StereographicProjection.html
 */
template<class T = double>
std::tuple<T,T> stereographic_proj(T twophi_crys, T twotheta_crys, T rad)
{
	const T sth = std::sin(twotheta_crys);
	const T cth = std::cos(twotheta_crys);
	const T sph = std::sin(twophi_crys);
	const T cph = std::cos(twophi_crys);

	T x = -T(2) * rad * sph * cth / (T(1) + cth*cph);
	T y = T(2) * rad * sth / (T(1) + cth*cph);

	return std::make_tuple(x, y);
}
// -----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// adapters
// ----------------------------------------------------------------------------
/**
 * vector-like access adapter to a matrix
 */
template<class t_mat> requires is_basic_mat<t_mat>
class matvec_adapter
{
public:
	using value_type = typename t_mat::value_type;
	using size_type = decltype(t_mat{}.size1());

public:
	matvec_adapter(const t_mat &mat) : m_mat{mat} {}
	~matvec_adapter() = default;

	size_type size() const { return m_mat.size1() * m_mat.size2(); }

	const value_type& operator[](size_type i) const
	{
		size_type row = i/m_mat.size2();
		size_type col = i%m_mat.size2();

		return m_mat(row, col);
	}

private:
	const t_mat& m_mat;
};


template<typename size_t, size_t N, typename T, template<size_t, size_t, class...> class t_mat_base>
class qvec_adapter : public t_mat_base<1, N, T>
{
public:
	// types
	using base_type = t_mat_base<1, N, T>;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qvec_adapter(const base_type& vec) : base_type{vec} {}

	static constexpr size_t size() { return N; }

	T& operator[](size_t i) { return base_type::operator()(i,0); }
	const T operator[](size_t i) const { return base_type::operator()(i,0); }
};

template<typename size_t, size_t ROWS, size_t COLS, typename T, template<size_t, size_t, class...> class t_mat_base>
class qmat_adapter : public t_mat_base<COLS, ROWS, T>
{
public:
	// types
	using base_type = t_mat_base<COLS, ROWS, T>;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qmat_adapter(const base_type& mat) : base_type{mat} {}

	static constexpr size_t size1() { return ROWS; }
	static constexpr size_t size2() { return COLS; }
};


template<typename size_t, size_t N, typename T, class t_vec_base>
class qvecN_adapter : public t_vec_base
{
public:
	// types
	using base_type = t_vec_base;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qvecN_adapter(const base_type& vec) : base_type{vec} {}

	static constexpr size_t size() { return N; }

	T& operator[](size_t i) { return static_cast<base_type&>(*this)[i]; }
	const T operator[](size_t i) const { return static_cast<const base_type&>(*this)[i]; }
};

template<typename size_t, size_t ROWS, size_t COLS, typename T, class t_mat_base>
class qmatNN_adapter : public t_mat_base
{
public:
	// types
	using base_type = t_mat_base;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qmatNN_adapter(const base_type& mat) : base_type{mat} {}

	// convert from a different matrix type
	template<class t_matOther> qmatNN_adapter(const t_matOther& matOther)
		requires is_basic_mat<t_matOther>
	{
		const std::size_t minRows = std::min(static_cast<std::size_t>(size1()), static_cast<std::size_t>(matOther.size1()));
		const std::size_t minCols = std::min(static_cast<std::size_t>(size2()), static_cast<std::size_t>(matOther.size2()));

		for(std::size_t i=0; i<minRows; ++i)
			for(std::size_t j=0; j<minCols; ++j)
				(*this)(i,j) = static_cast<value_type>(matOther(i,j));
	}

	static constexpr size_t size1() { return ROWS; }
	static constexpr size_t size2() { return COLS; }
};
// ----------------------------------------------------------------------------
}



namespace tl2_ops {
// ----------------------------------------------------------------------------
// vector operators
// ----------------------------------------------------------------------------

/**
 * unary +
 */
template<class t_vec>
const t_vec& operator+(const t_vec& vec1)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	return vec1;
}


/**
 * unary -
 */
template<class t_vec>
t_vec operator-(const t_vec& vec1)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	t_vec vec(vec1.size());

	for(std::size_t i=0; i<vec1.size(); ++i)
		vec[i] = -vec1[i];

	return vec;
}


/**
 * binary +
 */
template<class t_vec>
t_vec operator+(const t_vec& vec1, const t_vec& vec2)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	if constexpr(tl2::is_dyn_vec<t_vec>)
		assert((vec1.size() == vec2.size()));

	t_vec vec(vec1.size());

	for(std::size_t i=0; i<vec1.size(); ++i)
		vec[i] = vec1[i] + vec2[i];

	return vec;
}


/**
 * binary -
 */
template<class t_vec>
t_vec operator-(const t_vec& vec1, const t_vec& vec2)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	return vec1 + (-vec2);
}


/**
 * vector * scalar
 */
template<class t_vec>
t_vec operator*(const t_vec& vec1, typename t_vec::value_type d)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	t_vec vec(vec1.size());

	for(std::size_t i=0; i<vec1.size(); ++i)
		vec[i] = vec1[i] * d;

	return vec;
}


/**
 * scalar * vector
 */
template<class t_vec>
t_vec operator*(typename t_vec::value_type d, const t_vec& vec)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
	//&& !tl2::is_basic_mat<typename t_vec::value_type>	// hack!
{
	return vec * d;
}


/**
 * vector / scalar
 */
template<class t_vec>
t_vec operator/(const t_vec& vec1, typename t_vec::value_type d)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	t_vec vec(vec1.size());

	for(std::size_t i=0; i<vec1.size(); ++i)
		vec[i] = vec1[i] / d;

	return vec;
}


/**
 * vector += vector
 */
template<class t_vec>
t_vec& operator+=(t_vec& vec1, const t_vec& vec2)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	vec1 = vec1 + vec2;
	return vec1;
}


/**
 * vector -= vector
 */
template<class t_vec>
t_vec& operator-=(t_vec& vec1, const t_vec& vec2)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	vec1 = vec1 - vec2;
	return vec1;
}


/**
 * vector *= scalar
 */
template<class t_vec>
t_vec& operator*=(t_vec& vec1, typename t_vec::value_type d)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	vec1 = vec1 * d;
	return vec1;
}

/**
 * vector /= scalar
 */
template<class t_vec>
t_vec& operator/=(t_vec& vec1, typename t_vec::value_type d)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	vec1 = vec1 / d;
	return vec1;
}


/**
 * operator <<
 */
template<class t_vec>
std::ostream& operator<<(std::ostream& ostr, const t_vec& vec)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	const std::size_t N = vec.size();

	for(std::size_t i=0; i<N; ++i)
	{
		ostr << vec[i];
		if(i < N-1)
			ostr << TL2_COLSEP << " ";
	}

	return ostr;
}


/**
 * operator >>
 */
template<class t_vec>
std::istream& operator>>(std::istream& istr, t_vec& vec)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	vec.clear();

	std::string str;
	std::getline(istr, str);

	std::vector<std::string> vecstr;
	boost::split(vecstr, str,
		[](auto c) -> bool { return c==TL2_COLSEP; },
		boost::token_compress_on);

	for(auto& tok : vecstr)
	{
		boost::trim(tok);
		typename t_vec::value_type c = tl2::stoval<typename t_vec::value_type>(tok);
		vec.emplace_back(std::move(c));
	}

	return istr;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// matrix operators
// ----------------------------------------------------------------------------

/**
 * unary +
 */
template<class t_mat>
const t_mat& operator+(const t_mat& mat1)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	return mat1;
}


/**
 * unary -
 */
template<class t_mat>
t_mat operator-(const t_mat& mat1)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	t_mat mat(mat1.size1(), mat1.size2());

	for(std::size_t i=0; i<mat1.size1(); ++i)
		for(std::size_t j=0; j<mat1.size2(); ++j)
			mat(i,j) = -mat1(i,j);

	return mat;
}


/**
 * binary +
 */
template<class t_mat>
t_mat operator+(const t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat1.size1() == mat2.size1() && mat1.size2() == mat2.size2()));

	t_mat mat(mat1.size1(), mat1.size2());

	for(std::size_t i=0; i<mat1.size1(); ++i)
		for(std::size_t j=0; j<mat1.size2(); ++j)
			mat(i,j) = mat1(i,j) + mat2(i,j);

	return mat;
}


/**
 * binary -
 */
template<class t_mat>
t_mat operator-(const t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	return mat1 + (-mat2);
}


/**
 * matrix * scalar
 */
template<class t_mat>
t_mat operator*(const t_mat& mat1, typename t_mat::value_type d)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	t_mat mat(mat1.size1(), mat1.size2());

	for(std::size_t i=0; i<mat1.size1(); ++i)
		for(std::size_t j=0; j<mat1.size2(); ++j)
			mat(i,j) = mat1(i,j) * d;

	return mat;
}

/**
 * scalar * matrix
 */
template<class t_mat>
t_mat operator*(typename t_mat::value_type d, const t_mat& mat)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	return mat * d;
}


/**
 * matrix / scalar
 */
template<class t_mat>
t_mat operator/(const t_mat& mat, typename t_mat::value_type d)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	using T = typename t_mat::value_type;
	return mat * (T(1)/d);
}


/**
 * matrix-matrix product
 */
template<class t_mat>
t_mat operator*(const t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	return tl2::prod<t_mat>(mat1, mat2);
}


/**
 * matrix *= scalar
 */
template<class t_mat>
t_mat& operator*=(t_mat& mat1, typename t_mat::value_type d)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	mat1 = mat1 * d;
	return mat1;
}


/**
* matrix *= matrix
 */
template<class t_mat>
t_mat& operator*=(t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	mat1 = mat1 * mat2;
	return mat1;
}


/**
 * matrix += matrix
 */
template<class t_mat>
t_mat& operator+=(t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	mat1 = mat1 + mat2;
	return mat1;
}


/**
 * matrix -= matrix
 */
template<class t_mat>
t_mat& operator-=(t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	mat1 = mat1 - mat2;
	return mat1;
}


/**
 * matrix /= scalar
 */
template<class t_mat>
t_mat& operator/=(t_mat& mat1, typename t_mat::value_type d)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	mat1 = mat1 / d;
	return mat1;
}


/**
 * operator <<
 */
template<class t_mat>
std::ostream& operator<<(std::ostream& ostr, const t_mat& mat)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	const std::size_t ROWS = mat.size1();
	const std::size_t COLS = mat.size2();

	for(std::size_t row=0; row<ROWS; ++row)
	{
		for(std::size_t col=0; col<COLS; ++col)
		{
			ostr << mat(row, col);
			if(col < COLS-1)
				ostr << TL2_COLSEP << " ";
		}

		if(row < ROWS-1)
			ostr << TL2_ROWSEP << " ";
	}

	return ostr;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// mixed operators
// ----------------------------------------------------------------------------

/**
 * matrix-vector product: c_i = a_ij b_j
 */
template<class t_mat, class t_vec>
t_vec operator*(const t_mat& mat, const t_vec& vec)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
	&& tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size2() == vec.size()));
	else
		static_assert(t_mat::size2() == t_vec::size());


	t_vec vecRet(mat.size1());

	for(std::size_t row=0; row<mat.size1(); ++row)
	{
		vecRet[row] = typename t_vec::value_type{/*0*/};
		for(std::size_t col=0; col<mat.size2(); ++col)
		{
			auto elem = mat(row, col) * vec[col];
			vecRet[row] = vecRet[row] + elem;
		}
	}

	return vecRet;
}

// ----------------------------------------------------------------------------
}



namespace tl2 {
// ----------------------------------------------------------------------------
// vector and matrix containers
// ----------------------------------------------------------------------------

template<class T = double, template<class...> class t_cont = std::vector>
requires is_basic_vec<t_cont<T>> && is_dyn_vec<t_cont<T>>
class vec : public t_cont<T>
{
public:
	using value_type = T;
	using container_type = t_cont<T>;

	using container_type::container_type;
	using container_type::size;
	using container_type::operator[];

	using typename container_type::iterator;
	using typename container_type::const_iterator;
	using typename container_type::size_type;
	using typename container_type::difference_type;
	using typename container_type::allocator_type;

	~vec() = default;

	vec(const vec<T, t_cont>& other)
		: container_type{other}
	{}

	vec(vec<T, t_cont>&& other)
		: container_type{std::forward<vec<T, t_cont>&&>(other)}
	{}

	template<class T_other, template<class...> class t_cont_other>
	vec(const vec<T_other, t_cont_other>& other)
	{
		this->operator=<T_other, t_cont_other>(other);
	}

	vec<T, t_cont>& operator=(const vec<T, t_cont>& other)
	{
		*static_cast<container_type*>(this) = other;
		return *this;
	}

	const vec<T, t_cont>& operator=(const vec<T, t_cont>& other) const
	{
		*static_cast<container_type*>(this) = other;
		return *this;
	}

	template<class T_other, template<class...> class t_cont_other>
	vec<T, t_cont>& operator=(const vec<T_other, t_cont_other>& other)
	{
		*this = convert<vec<T, t_cont>, vec<T_other, t_cont_other>>(other);
		return *this;
	}

	template<class T_other, template<class...> class t_cont_other>
	const vec<T, t_cont>& operator=(const vec<T_other, t_cont_other>& other) const
	{
		*this = convert<vec<T, t_cont>, vec<T_other, t_cont_other>>(other);
		return *this;
	}

	vec(std::size_t SIZE, const T* arr = nullptr) : container_type(SIZE)
	{
		if(arr)
			from_array(arr);
	}


	const value_type& operator()(std::size_t i) const { return this->operator[](i); }
	value_type& operator()(std::size_t i) { return this->operator[](i); }


	void from_array(const T* arr)
	{
		// initialise from given array data
		for(std::size_t i=0; i<size(); ++i)
			this->operator[](i) = arr[i];
	}

	void to_array(T* arr) const
	{
		// write elements to array
		for(std::size_t i=0; i<size(); ++i)
			arr[i] = this->operator[](i);
	}

	friend vec operator+(const vec& vec1, const vec& vec2) { return tl2_ops::operator+(vec1, vec2); }
	friend vec operator-(const vec& vec1, const vec& vec2) { return tl2_ops::operator-(vec1, vec2); }
	friend const vec& operator+(const vec& vec1) { return tl2_ops::operator+(vec1); }
	friend vec operator-(const vec& vec1) { return tl2_ops::operator-(vec1); }

	friend vec operator*(value_type d, const vec& vec1) { return tl2_ops::operator*(d, vec1); }
	friend vec operator*(const vec& vec1, value_type d) { return tl2_ops::operator*(vec1, d); }
	friend vec operator/(const vec& vec1, value_type d) { return tl2_ops::operator/(vec1, d); }

	vec& operator*=(const vec& vec2) { return tl2_ops::operator*=(*this, vec2); }
	vec& operator+=(const vec& vec2) { return tl2_ops::operator+=(*this, vec2); }
	vec& operator-=(const vec& vec2) { return tl2_ops::operator-=(*this, vec2); }
	vec& operator*=(value_type d) { return tl2_ops::operator*=(*this, d); }
	vec& operator/=(value_type d) { return tl2_ops::operator/=(*this, d); }
};


template<class T=double, template<class...> class t_cont = std::vector>
requires is_basic_vec<t_cont<T>> && is_dyn_vec<t_cont<T>>
class mat
{
public:
	using value_type = T;
	using container_type = t_cont<T>;

	mat() = default;
	~mat() = default;

	mat(std::size_t ROWS, std::size_t COLS, const T* arr = nullptr)
		: m_data(ROWS*COLS), m_rowsize{ROWS}, m_colsize{COLS}
	{
		if(arr)
			from_array(arr);
	}


	template<class T_other, template<class...> class t_cont_other>
	mat(const mat<T_other, t_cont_other>& other)
	{
		this->operator=<T_other, t_cont_other>(other);
	}

	template<class T_other, template<class...> class t_cont_other>
	mat<T, t_cont>& operator=(const vec<T_other, t_cont_other>& other)
	{
		*this = convert<mat<T, t_cont>, mat<T_other, t_cont_other>>(other);

		this->m_rowsize = other.m_rowsize;
		this->m_colsize = other.m_colsize;

		return *this;
	}

	template<class T_other, template<class...> class t_cont_other>
	const mat<T, t_cont>& operator=(const vec<T_other, t_cont_other>& other) const
	{
		*this = convert<mat<T, t_cont>, mat<T_other, t_cont_other>>(other);

		this->m_rowsize = other.m_rowsize;
		this->m_colsize = other.m_colsize;

		return *this;
	}


	std::size_t size1() const { return m_rowsize; }
	std::size_t size2() const { return m_colsize; }


	// element access
	const T& operator()(std::size_t row, std::size_t col) const { return m_data[row*m_colsize + col]; }
	T& operator()(std::size_t row, std::size_t col) { return m_data[row*m_colsize + col]; }


	void from_array(const T* arr)
	{
		// initialise from given array data
		for(std::size_t i=0; i<m_rowsize; ++i)
			for(std::size_t j=0; j<m_colsize; ++j)
				this->operator()(i,j) = arr[i*m_colsize + j];
	}

	void to_array(T* arr) const
	{
		// write elements to array
		for(std::size_t i=0; i<m_rowsize; ++i)
			for(std::size_t j=0; j<m_colsize; ++j)
				arr[i*m_colsize + j] = this->operator()(i,j);
	}

	friend mat operator+(const mat& mat1, const mat& mat2) { return tl2_ops::operator+(mat1, mat2); }
	friend mat operator-(const mat& mat1, const mat& mat2) { return tl2_ops::operator-(mat1, mat2); }
	friend const mat& operator+(const mat& mat1) { return tl2_ops::operator+(mat1); }
	friend mat operator-(const mat& mat1) { return tl2_ops::operator-(mat1); }

	friend mat operator*(const mat& mat1, const mat& mat2) { return tl2_ops::operator*(mat1, mat2); }
	friend mat operator*(const mat& mat1, value_type d) { return tl2_ops::operator*(mat1, d); }
	friend mat operator*(value_type d, const mat& mat1) { return tl2_ops::operator*(d, mat1); }
	friend mat operator/(const mat& mat1, value_type d) { return tl2_ops::operator/(mat1, d); }

	template<class t_vec> requires is_basic_vec<t_cont<T>> && is_dyn_vec<t_cont<T>>
	friend t_vec operator*(const mat& mat1, const t_vec& vec2) { return tl2_ops::operator*(mat1, vec2); }

	mat& operator*=(const mat& mat2) { return tl2_ops::operator*=(*this, mat2); }
	mat& operator+=(const mat& mat2) { return tl2_ops::operator+=(*this, mat2); }
	mat& operator-=(const mat& mat2) { return tl2_ops::operator-=(*this, mat2); }
	mat& operator*=(value_type d) { return tl2_ops::operator*=(*this, d); }
	mat& operator/=(value_type d) { return tl2_ops::operator/=(*this, d); }

private:
	container_type m_data{};
	std::size_t m_rowsize{0};
	std::size_t m_colsize{0};
};

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// scalar algos
// ----------------------------------------------------------------------------

/**
 * are two scalars equal within an epsilon range?
 */
template<class T>
bool equals(T t1, T t2, T eps = std::numeric_limits<T>::epsilon())
requires is_scalar<T>
{
	return std::abs(t1 - t2) <= eps;
}


/**
 * is the given value an integer?
 */
template<class T = double>
bool is_integer(T val, T eps = std::numeric_limits<T>::epsilon())
requires is_scalar<T>
{
	T val_diff = val - std::round(val);
	return equals<T>(val_diff, T(0), eps);
}


/**
 * get next multiple of the given granularity
 */
template<typename t_num = unsigned int>
t_num next_multiple(t_num num, t_num granularity)
requires is_scalar<t_num>
{
	t_num div = num / granularity;
	bool rest_is_0 = 1;

	if constexpr(std::is_floating_point_v<t_num>)
	{
		div = std::floor(div);
		t_num rest = std::fmod(num, granularity);
		rest_is_0 = equals(rest, t_num{0});
	}
	else
	{
		t_num rest = num % granularity;
		rest_is_0 = (rest==0);
	}

	return rest_is_0 ? num : (div+1) * granularity;
}


/**
 * mod operation, keeping result positive
 */
template<class t_real>
t_real mod_pos(t_real val, t_real tomod = t_real{2}*pi<t_real>)
requires is_scalar<t_real>
{
	val = std::fmod(val, tomod);
	if(val < t_real(0))
		val += tomod;

	return val;
}

/**
 * are two angles equal within an epsilon range?
 */
template<class T>
bool angle_equals(T t1, T t2,
	T eps = std::numeric_limits<T>::epsilon(),
	T tomod = T{2}*pi<T>)
requires is_scalar<T>
{
	t1 = mod_pos<T>(t1, tomod);
	t2 = mod_pos<T>(t2, tomod);

	return std::abs(t1 - t2) <= eps;
}


/**
 * are two complex numbers equal within an epsilon range?
 */
template<class T> requires is_complex<T>
bool equals(const T& t1, const T& t2,
	typename T::value_type eps = std::numeric_limits<typename T::value_type>::epsilon())
{
	return (std::abs(t1.real() - t2.real()) <= eps) &&
		(std::abs(t1.imag() - t2.imag()) <= eps);
}


/**
 * are two complex numbers equal within an epsilon range?
 */
template<class T> requires is_complex<T>
bool equals(const T& t1, const T& t2, const T& eps)
{
	return equals<T>(t1, t2, eps.real());
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// n-dim algos
// ----------------------------------------------------------------------------

/**
 * are two vectors equal within an epsilon range?
 */
template<class t_vec, class t_num = typename t_vec::value_type>
bool equals(const t_vec& vec1, const t_vec& vec2,
	t_num eps = std::numeric_limits<t_num>::epsilon(),
	int _maxSize = -1)
requires is_basic_vec<t_vec>
{
	// size has to be equal
	if(vec1.size() != vec2.size())
		return false;

	std::size_t maxSize = vec1.size();
	if(_maxSize >= 0)
		maxSize = std::min(std::size_t(_maxSize), maxSize);

	// check each element
	for(std::size_t i=0; i<maxSize; ++i)
	{
		if constexpr(tl2::is_complex<t_num>)
		{
			if(!equals<t_num>(vec1[i], vec2[i], eps.real()))
				return false;
		}
		else if constexpr(tl2::is_scalar<t_num>)
		{
			if(!equals<t_num>(vec1[i], vec2[i], eps))
				return false;
		}
	}

	return true;
}


/**
 * are two matrices equal within an epsilon range?
 */
template<class t_mat, class t_real>
bool equals(const t_mat& mat1, const t_mat& mat2,
	t_real eps = std::numeric_limits<t_real>::epsilon(),
	int _maxSize = -1)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;

	if(mat1.size1() != mat2.size1() || mat1.size2() != mat2.size2())
		return false;

	std::size_t maxSize1 = mat1.size1();
	std::size_t maxSize2 = mat1.size2();
	if(_maxSize >= 0)
	{
		maxSize1 = std::min(std::size_t(_maxSize), maxSize1);
		maxSize2 = std::min(std::size_t(_maxSize), maxSize2);
	}

	for(std::size_t i=0; i<maxSize1; ++i)
	{
		for(std::size_t j=0; j<maxSize2; ++j)
		{
			if constexpr(is_complex<decltype(eps)>)
			{
				if(!equals<T>(mat1(i,j), mat2(i,j), eps.real()))
					return false;
			}
			else
			{
				if(!equals<T>(mat1(i,j), mat2(i,j), eps))
					return false;
			}
		}
	}

	return true;
}


/**
 * check if two collections of matrices or vectors are equal
 */
template<class t_obj, template<class...> class t_vec = std::vector>
bool equals_all(const t_vec<t_obj>& vec1, const t_vec<t_obj>& _vec2,
	typename t_obj::value_type eps = std::numeric_limits<typename t_obj::value_type>::epsilon(),
	int maxSize=-1)
{
	auto vec2 = _vec2;
	if(vec1.size() != vec2.size())
		return false;

	for(const auto& obj1 : vec1)
	{
		// find obj1 in vec2
		auto iter = std::find_if(vec2.crbegin(), vec2.crend(),
		[&obj1, eps, maxSize](const t_obj& obj2) -> bool
		{
			return tl2::equals<t_obj>(obj1, obj2, eps, maxSize);
		});

		// not found
		if(iter == vec2.crend())
			return false;

		// remove already checked element
		vec2.erase(iter.base()-1);
	}

	return true;
}


/**
 * remove duplicate vectors or matrices in the container
 */
template<class t_obj, template<class...> class t_cont = std::vector>
t_cont<t_obj> remove_duplicates(const t_cont<t_obj>& objs,
	typename t_obj::value_type eps = std::numeric_limits<typename t_obj::value_type>::epsilon())
{
	t_cont<t_obj> newobjs;
	newobjs.reserve(objs.size());

	for(const auto& elem : objs)
	{
		// find obj in container
		auto iter = std::find_if(newobjs.cbegin(), newobjs.cend(),
		[&elem, eps](const t_obj& elem2) -> bool
		{
			return tl2::equals<t_obj>(elem, elem2, eps);
		});

		// not found
		if(iter == newobjs.cend())
			newobjs.push_back(elem);
	}

	return newobjs;
}


/**
 * set submatrix to unit
 */
template<class t_mat>
void unit(t_mat& mat, std::size_t rows_begin, std::size_t cols_begin, std::size_t rows_end, std::size_t cols_end)
requires is_basic_mat<t_mat>
{
	for(std::size_t i=rows_begin; i<rows_end; ++i)
		for(std::size_t j=cols_begin; j<cols_end; ++j)
			mat(i,j) = (i==j ? 1 : 0);
}


/**
 * unit matrix
 */
template<class t_mat>
t_mat unit(std::size_t N1, std::size_t N2)
requires is_basic_mat<t_mat>
{
	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	unit<t_mat>(mat, 0,0, mat.size1(),mat.size2());
	return mat;
}


/**
 * unit matrix
 */
template<class t_mat>
t_mat unit(std::size_t N=0)
requires is_basic_mat<t_mat>
{
	return unit<t_mat>(N,N);
}


/**
 * is mat a unit matrix?
 */
template<class t_mat, class t_scalar = typename t_mat::value_type>
bool is_unit(const t_mat& mat,
	t_scalar eps = std::numeric_limits<t_scalar>::epsilon())
requires is_mat<t_mat>
{
	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			if(i==j && !equals<t_scalar>(mat(i, j), t_scalar(1), eps))
				return false;
			if(i!=j && !equals<t_scalar>(mat(i, j), t_scalar(0), eps))
				return false;
		}
	}

	return true;
}


/**
 * zero matrix
 */
template<class t_mat>
t_mat zero(std::size_t N1, std::size_t N2)
requires is_basic_mat<t_mat>
{
	using size_t = decltype(t_mat{}.size1());

	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	for(size_t i=0; i<mat.size1(); ++i)
		for(size_t j=0; j<mat.size2(); ++j)
			mat(i,j) = 0;

	return mat;
}


/**
 * zero matrix
 */
template<class t_mat>
t_mat zero(std::size_t N=0)
requires is_basic_mat<t_mat>
{
	return zero<t_mat>(N, N);
}


/**
 * is mat a zero matrix?
 */
template<class t_mat, class t_scalar = typename t_mat::value_type>
bool is_zero(const t_mat& mat,
	t_scalar eps = std::numeric_limits<t_scalar>::epsilon())
requires is_mat<t_mat>
{
	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			if(!equals<t_scalar>(mat(i, j), t_scalar(0), eps))
				return false;
		}
	}

	return true;
}


/**
 * zero vector
 */
template<class t_vec>
t_vec zero(std::size_t N /* = 0*/)
requires is_basic_vec<t_vec>
{
	using size_t = decltype(t_vec{}.size());

	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(N);

	for(size_t i=0; i<vec.size(); ++i)
		vec[i] = 0;

	return vec;
}


/**
 * is vec a zero vector?
 */
template<class t_vec, class t_scalar = typename t_vec::value_type>
bool is_zero(const t_vec& vec,
	t_scalar eps = std::numeric_limits<t_scalar>::epsilon())
requires is_vec<t_vec>
{
	for(std::size_t i=0; i<vec.size(); ++i)
	{
		if(!equals<t_scalar>(vec[i], t_scalar(0), eps))
			return false;
	}

	return true;
}


/**
 * random scalar
 */
template<typename t_scalar>
t_scalar rand()
requires is_scalar<t_scalar> && (!is_complex<t_scalar>)
{
	return get_rand<t_scalar>(0, 1);
};


/**
 * random complex number
 */
template<typename t_cplx, class t_scalar = typename t_cplx::value_type>
t_cplx rand()
requires is_complex<t_cplx> && is_scalar<t_scalar>
{
	return t_cplx(rand<t_scalar>(), rand<t_scalar>());
};


/**
 * random vector
 */
template<class t_vec>
t_vec rand(std::size_t N)
requires is_basic_vec<t_vec>
{
	using size_t = decltype(t_vec{}.size());
	using t_scalar = typename t_vec::value_type;

	t_vec vec{};
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(N);

	for(size_t i=0; i<vec.size(); ++i)
		vec[i] = rand<t_scalar>();

	return vec;
}


/**
 * random matrix
 */
template<class t_mat>
t_mat rand(std::size_t N1, std::size_t N2)
requires is_basic_mat<t_mat>
{
	using size_t = decltype(t_mat{}.size1());
	using t_scalar = typename t_mat::value_type;

	t_mat mat{};
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	for(size_t i=0; i<mat.size1(); ++i)
		for(size_t j=0; j<mat.size2(); ++j)
			mat(i,j) = rand<t_scalar>();

	return mat;
}


/**
 * row permutation matrix
 * @see https://en.wikipedia.org/wiki/Permutation_matrix
 */
template<class t_mat>
t_mat perm(std::size_t N1, std::size_t N2, std::size_t from, std::size_t to)
requires is_basic_mat<t_mat>
{
	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	unit<t_mat>(mat, 0,0, mat.size1(),mat.size2());

	mat(from, from) = mat(to, to) = 0;
	mat(from, to) = mat(to, from) = 1;

	return mat;
}


/**
 * create matrix representation of a permutation
 * @see https://en.wikipedia.org/wiki/Permutation_matrix
 */
template<class t_mat, class t_perm=std::vector<std::size_t>>
t_mat perm(const t_perm& perm)
requires is_basic_mat<t_mat>
{
	auto N = perm.size();
	t_mat mat = zero<t_mat>(N, N);

	for(decltype(N) i=0; i<perm.size(); ++i)
		mat(i, perm[i]) = mat(perm[i], i) = 1;

	return mat;
}


/**
 * diagonal matrix
 */
template<class t_mat, class t_vec = std::vector<typename t_mat::value_type>>
t_mat diag(const t_vec& vals)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	const std::size_t N = vals.size();
	t_mat mat = zero<t_mat>(N);

	// static matrix does not necessarily have the required size!
	if constexpr(!tl2::is_dyn_mat<t_mat>)
		assert(mat.size1() == mat.size1() && mat.size1() == N);

	for(std::size_t i=0; i<std::min(mat.size1(), N); ++i)
		mat(i,i) = vals[i];

	return mat;
}


/**
 * diagonal matrix
 */
template<class t_mat, class t_val = typename t_mat::value_type>
t_mat diag(std::size_t N, const t_val& val)
requires is_basic_mat<t_mat>
{
	t_mat mat = zero<t_mat>(N, N);

	for(std::size_t i=0; i<N; ++i)
		mat(i, i) = val;

	return mat;
}


/**
 * vector of diagonal matrix elements
 */
template<class t_vec, class t_mat>
t_vec diag_vec(const t_mat& mat)
requires is_vec<t_vec> && is_mat<t_mat>
{
	std::size_t N = std::min(mat.size1(), mat.size2());

	t_vec vec = zero<t_vec>(N);
	for(std::size_t i=0; i<N; ++i)
		vec[i] = mat(i,i);

	return vec;
}


/**
 * tests for zero scalar
 */
template<class t_real>
bool equals_0(t_real val, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_scalar<t_real>
{
	return equals<t_real>(val, t_real(0), eps);
}


/**
 * tests for zero vector
 */
template<class t_vec>
bool equals_0(const t_vec& vec,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
requires is_basic_vec<t_vec>
{
	return equals<t_vec>(vec, zero<t_vec>(vec.size()), eps);
}


/**
 * tests for zero matrix
 */
template<class t_mat>
bool equals_0(const t_mat& mat,
	typename t_mat::value_type eps = std::numeric_limits<typename t_mat::value_type>::epsilon())
requires is_mat<t_mat>
{
	return equals<t_mat>(mat, zero<t_mat>(mat.size1(), mat.size2()), eps);
}


/**
 * tests for symmetric or hermitian matrix
 */
template<class t_mat, class t_real = typename t_mat::value_type>
bool is_symm_or_herm(const t_mat& mat,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat>
{
	using t_elem = typename t_mat::value_type;
	if(mat.size1() != mat.size2())
		return false;

	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=i+1; j<mat.size2(); ++j)
		{
			if constexpr(is_complex<t_elem>)
			{
				// not hermitian?
				if(!equals<t_elem>(mat(i,j), std::conj(mat(j,i)), eps))
					return false;
			}
			else
			{
				// not symmetric?
				if(!equals<t_elem>(mat(i,j), mat(j,i), eps))
					return false;
			}
		}
	}

	return true;
}


/**
 * tests for skew-symmetric or skew-hermitian matrix
 * @see https://en.wikipedia.org/wiki/Skew-Hermitian_matrix
 */
template<class t_mat, class t_real = typename t_mat::value_type>
bool is_skew_symm_or_herm(const t_mat& mat,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat>
{
	using t_elem = typename t_mat::value_type;
	if(mat.size1() != mat.size2())
		return false;

	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=i+1; j<mat.size2(); ++j)
		{
			if constexpr(is_complex<t_elem>)
			{
				// not hermitian?
				if(!equals<t_elem>(mat(i,j), -std::conj(mat(j,i)), eps))
					return false;
			}
			else
			{
				// not skew-symmetric?
				if(!equals<t_elem>(mat(i,j), -mat(j,i), eps))
					return false;
			}
		}
	}

	return true;
}


/**
 * tests for diagonal matrix
 */
template<class t_mat, class t_real = typename t_mat::value_type>
bool is_diag(const t_mat& mat,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat>
{
	using t_elem = typename t_mat::value_type;
	if(mat.size1() != mat.size2())
		return false;

	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			if(i == j)
				continue;

			if constexpr(is_complex<t_elem>)
			{
				if(!equals<t_elem>(mat(i,j), t_elem(0, 0), eps))
					return false;
			}
			else
			{
				if(!equals<t_elem>(mat(i,j), t_elem(0.), eps))
					return false;
			}
		}
	}

	return true;
}


/**
 * transpose matrix
 * WARNING: not possible for static non-square matrix!
 */
template<class t_mat>
t_mat trans(const t_mat& mat)
requires is_mat<t_mat>
{
	using t_idxtype = decltype(mat.size1());

	t_mat mat2;
	if constexpr(is_dyn_mat<t_mat>)
		mat2 = t_mat(mat.size2(), mat.size1());

	for(t_idxtype i=0; i<mat.size1(); ++i)
		for(t_idxtype j=0; j<mat.size2(); ++j)
			mat2(j,i) = mat(i,j);

	return mat2;
}


// -----------------------------------------------------------------------------
/**
 * set values lower than epsilon to zero
 * scalar version
 */
template<typename t_real>
void set_eps_0(t_real& d, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_scalar<t_real> && (!is_complex<t_real>)
{
	if(std::abs(d) < eps)
		d = t_real(0);
};


template<typename t_real, typename t_real_2>
void set_eps_0(t_real& d, t_real_2 eps = std::numeric_limits<t_real_2>::epsilon())
requires is_scalar<t_real> && is_scalar<t_real_2> && (!is_complex<t_real> && !is_complex<t_real_2>)
{
	if(std::abs(d) < t_real(eps))
		d = t_real(0);
};


/**
 * set values lower than epsilon to zero
 * complex version
 */
template<typename t_cplx, class t_real = typename t_cplx::value_type>
void set_eps_0(t_cplx& c, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_complex<t_cplx> && is_scalar<t_real>
{
	if(std::abs(c.real()) < eps)
		c.real(t_real(0));
	if(std::abs(c.imag()) < eps)
		c.imag(t_real(0));
};


/**
 * set values lower than epsilon to zero
 * vector version
 */
template<typename t_vec, typename t_real = typename t_vec::value_type>
void set_eps_0(t_vec& vec, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_basic_vec<t_vec>
{
	using t_elem = typename t_vec::value_type;

	for(std::size_t i=0; i<vec.size(); ++i)
		set_eps_0<t_elem>(vec[i], eps);
};


/**
 * set values lower than epsilon to zero
 * matrix version
 */
template<typename t_mat, typename t_real = typename t_mat::value_type>
void set_eps_0(t_mat& mat, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_basic_mat<t_mat>
{
	using t_elem = typename t_mat::value_type;

	for(std::size_t i=0; i<mat.size1(); ++i)
		for(std::size_t j=0; j<mat.size2(); ++j)
			set_eps_0<t_elem>(mat(i,j), eps);
};
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
/**
 * prints a matrix in a 'nicely' formatted form
 */
template<class t_mat, class t_elem = typename t_mat::value_type, class t_real = double>
std::ostream& niceprint(std::ostream& ostr, const t_mat& mat,
	t_real eps = 1e-6, unsigned int prec = 6)
requires tl2::is_basic_mat<t_mat> /*&& tl2::is_dyn_mat<t_mat>*/
{
	ostr.precision(prec);

	const std::size_t ROWS = mat.size1();
	const std::size_t COLS = mat.size2();

	for(std::size_t i=0; i<ROWS; ++i)
	{
		ostr << "(";
		for(std::size_t j=0; j<COLS; ++j)
		{
			t_elem elem = mat(i, j);
			set_eps_0(elem, eps);

			unsigned int field_width = prec * 1.5 + 4;
			if constexpr(is_complex<t_elem>)
				field_width *= 3;
			ostr << std::setw(field_width) << std::right << elem;
		}
		ostr << ")\n";
	}

	return ostr;
}


/**
 * prints a vector in a 'nicely' formatted form
 */
template<class t_vec, class t_elem = typename t_vec::value_type, class t_real = double>
std::ostream& niceprint(std::ostream& ostr, const t_vec& vec,
	t_real eps = 1e-6, unsigned int prec = 6)
requires tl2::is_basic_vec<t_vec> /*&& tl2::is_dyn_vec<t_vec>*/
{
	ostr.precision(prec);

	const std::size_t ROWS = vec.size();

	for(std::size_t i=0; i<ROWS; ++i)
	{
		ostr << "(";
		t_elem elem = vec[i];
		set_eps_0(elem, eps);

		unsigned int field_width = prec * 1.5 + 4;
		if constexpr(is_complex<t_elem>)
			field_width *= 3;
		ostr << std::setw(field_width) << std::right << elem;
		ostr << ")\n";
	}

	return ostr;
}
// -----------------------------------------------------------------------------



/**
 * create a vector with given size if it is dynamic
 */
template<class t_vec>
t_vec create(std::size_t size=3)
requires is_basic_vec<t_vec>
{
	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(size);

	return vec;
}


/**
 * create a matrix with given sizes if it is dynamic
 */
template<class t_mat>
t_mat create(std::size_t size1, std::size_t size2)
requires is_basic_mat<t_mat>
{
	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat{size1, size2};

	return mat;
}


/**
 * create vector from initializer_list
 */
template<class t_vec, template<class...> class t_cont = std::initializer_list>
t_vec create(const t_cont<typename t_vec::value_type>& lst)
requires is_basic_vec<t_vec>
{
	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(lst.size());

	auto iterLst = lst.begin();
	auto size = vec.size();
	using local_size_t = std::decay_t<decltype(size)>;
	for(local_size_t i=0; i<size; ++i)
	{
		if(iterLst != lst.end())
		{
			vec[i] = *iterLst;
			std::advance(iterLst, 1);
		}
		else	// vector larger than given list?
		{
			vec[i] = 0;
		}
	}

	return vec;
}


/**
 * create matrix from nested initializer_lists in columns/rows order
 */
template<class t_mat,
	template<class...> class t_cont_outer = std::initializer_list,
	template<class...> class t_cont = std::initializer_list>
t_mat create_mat(const t_cont_outer<t_cont<typename t_mat::value_type>>& lst)
requires is_mat<t_mat>
{
	const std::size_t iCols = lst.size();
	const std::size_t iRows = lst.begin()->size();

	t_mat mat = unit<t_mat>(iRows, iCols);

	auto iterCol = lst.begin();
	for(std::size_t iCol=0; iCol<iCols; ++iCol)
	{
		auto iterRow = iterCol->begin();
		for(std::size_t iRow=0; iRow<iRows; ++iRow)
		{
			mat(iRow, iCol) = *iterRow;
			std::advance(iterRow, 1);
		}

		std::advance(iterCol, 1);
	}

	return mat;
}


/**
 * create matrix from nested initializer_lists in columns/rows order
 */
template<class t_mat,
	template<class...> class t_cont_outer = std::initializer_list,
	template<class...> class t_cont = std::initializer_list>
t_mat create(const t_cont_outer<t_cont<typename t_mat::value_type>>& lst)
requires is_mat<t_mat>
{
	return create_mat<t_mat, t_cont_outer, t_cont>(lst);
}


/**
 * create matrix from column (or row) vectors
 */
template<class t_mat, class t_vec, template<class...> class t_cont_outer = std::initializer_list>
t_mat create(const t_cont_outer<t_vec>& lst, bool bRow = false)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	const std::size_t iCols = lst.size();
	const std::size_t iRows = lst.begin()->size();

	t_mat mat = unit<t_mat>(iRows, iCols);

	auto iterCol = lst.begin();
	for(std::size_t iCol=0; iCol<iCols; ++iCol)
	{
		for(std::size_t iRow=0; iRow<iRows; ++iRow)
			mat(iRow, iCol) = (*iterCol)[iRow];
		std::advance(iterCol, 1);
	}

	if(bRow) mat = trans<t_mat>(mat);
	return mat;
}


/**
 * create matrix from initializer_list in column/row order
 */
template<class t_mat>
t_mat create_mat(const std::initializer_list<typename t_mat::value_type>& lst)
requires is_mat<t_mat>
{
	const std::size_t N = std::sqrt(lst.size());

	t_mat mat = unit<t_mat>(N, N);

	auto iter = lst.begin();
	for(std::size_t iRow=0; iRow<N; ++iRow)
	{
		for(std::size_t iCol=0; iCol<N; ++iCol)
		{
			mat(iRow, iCol) = *iter;
			std::advance(iter, 1);
		}
	}

	return mat;
}


/**
 * create matrix from initializer_list in column/row order
 */
template<class t_mat>
t_mat create(const std::initializer_list<typename t_mat::value_type>& lst)
requires is_mat<t_mat>
{
	return create_mat<t_mat>(lst);
}


/**
 * linearise a matrix to a vector container
 */
template<class t_vec, class t_mat>
t_vec convert(const t_mat& mat)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	using T_dst = typename t_vec::value_type;
	using t_idx = decltype(mat.size1());

	t_vec vec;
	vec.reserve(mat.size1()*mat.size2());

	for(t_idx iRow=0; iRow<mat.size1(); ++iRow)
		for(t_idx iCol=0; iCol<mat.size2(); ++iCol)
			vec.push_back(T_dst(mat(iRow, iCol)));

	return vec;
}


/**
 * converts matrix containers of different value types
 */
template<class t_mat_dst, class t_mat_src>
t_mat_dst convert(const t_mat_src& mat)
requires is_mat<t_mat_dst> && is_mat<t_mat_src>
{
	using T_dst = typename t_mat_dst::value_type;
	using t_idx = decltype(mat.size1());

	// in case the static size of the destination vector is larger than the source's
	t_idx maxRows = std::max(mat.size1(), t_idx(t_mat_dst{}.size1()));
	t_idx maxCols = std::max(mat.size2(), t_idx(t_mat_dst{}.size2()));

	t_mat_dst matdst = unit<t_mat_dst>(maxRows, maxCols);

	for(t_idx iRow=0; iRow<mat.size1(); ++iRow)
		for(t_idx iCol=0; iCol<mat.size2(); ++iCol)
			matdst(iRow, iCol) = T_dst(mat(iRow, iCol));

	return matdst;
}


/**
 * converts vector containers of different value types
 */
template<class t_vec_dst, class t_vec_src>
t_vec_dst convert(const t_vec_src& vec)
requires is_vec<t_vec_dst> && is_vec<t_vec_src>
{
	using T_dst = typename t_vec_dst::value_type;
	using t_idx = decltype(vec.size());

	t_vec_dst vecdst = create<t_vec_dst>(vec.size());

	for(t_idx i=0; i<vec.size(); ++i)
		vecdst[i] = T_dst(vec[i]);

	return vecdst;
}


/**
 * converts a container of objects
 */
template<class t_obj_dst, class t_obj_src, template<class...> class t_cont>
t_cont<t_obj_dst> convert(const t_cont<t_obj_src>& src_objs)
requires (is_vec<t_obj_dst> || is_mat<t_obj_dst>) && (is_vec<t_obj_src> || is_mat<t_obj_src>)
{
	t_cont<t_obj_dst> dst_objs;
	dst_objs.reserve(src_objs.size());

	for(const t_obj_src& src_obj : src_objs)
		dst_objs.emplace_back(convert<t_obj_dst, t_obj_src>(src_obj));

	return dst_objs;
}


/**
 * get a column vector from a matrix
 */
template<class t_mat, class t_vec>
t_vec col(const t_mat& mat, std::size_t col)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(mat.size1());

	for(std::size_t i=0; i<mat.size1(); ++i)
		vec[i] = mat(i, col);

	return vec;
}


/**
 * get a row vector from a matrix
 */
template<class t_mat, class t_vec>
t_vec row(const t_mat& mat, std::size_t row)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(mat.size2());

	auto size = std::min(mat.size2(), vec.size());
	for(std::size_t i=0; i<std::size_t(size); ++i)
		vec[i] = mat(row, i);

	return vec;
}


/**
 * set a column vector in a matrix
 */
template<class t_mat, class t_vec>
void set_col(t_mat& mat, const t_vec& vec, std::size_t col)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	for(std::size_t i=0; i<std::min<std::size_t>(mat.size1(), vec.size()); ++i)
		mat(i, col) = vec[i];
}


/**
 * set a row vector in a matrix
 */
template<class t_mat, class t_vec>
void set_row(t_mat& mat, const t_vec& vec, std::size_t row)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	for(std::size_t i=0; i<std::min<std::size_t>(mat.size1(), vec.size()); ++i)
		mat(row, i) = vec[i];
}


/**
 * inner product <vec1|vec2>
 */
template<class t_vec>
typename t_vec::value_type inner(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec>
{
	typename t_vec::value_type val{0};
	auto size = vec1.size();
	using local_size_t = std::decay_t<decltype(size)>;

	for(local_size_t i=0; i<size; ++i)
	{
		if constexpr(is_complex<typename t_vec::value_type>)
			val += std::conj(vec1[i]) * vec2[i];
		else
			val += vec1[i] * vec2[i];
	}

	return val;
}


/**
 * inner product <vec1|vec2> (without conjugation of complex vector)
 */
template<class t_vec>
typename t_vec::value_type inner_noconj(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec>
{
	typename t_vec::value_type val{0};
	auto size = vec1.size();
	using local_size_t = std::decay_t<decltype(size)>;

	for(local_size_t i=0; i<size; ++i)
		val += vec1[i] * vec2[i];

	return val;
}


/**
 * inner product between two vectors of different type
 */
template<class t_vec1, class t_vec2>
typename t_vec1::value_type inner(const t_vec1& vec1, const t_vec2& vec2)
requires is_basic_vec<t_vec1> && is_basic_vec<t_vec2>
{
	if(vec1.size()==0 || vec2.size()==0)
		return typename t_vec1::value_type{};

	// first element
	auto val = vec1[0]*vec2[0];

	// remaining elements
	for(std::size_t i=1; i<std::min(vec1.size(), vec2.size()); ++i)
	{
		if constexpr(is_complex<typename t_vec1::value_type>)
		{
			auto prod = std::conj(vec1[i]) * vec2[i];
			val = val + prod;
		}
		else
		{
			auto prod = vec1[i]*vec2[i];
			val = val + prod;
		}
	}

	return val;
}


/**
 * matrix-matrix product: c_ij = a_ik b_kj
 */
template<class t_mat>
t_mat prod(const t_mat& mat1, const t_mat& mat2, bool assert_sizes/*=true*/)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	// if not asserting sizes, the inner size will use the minimum of the two matrix sizes
	if(assert_sizes)
	{
		if constexpr(tl2::is_dyn_mat<t_mat>)
			assert((mat1.size2() == mat2.size1()));
		else
			static_assert(t_mat::size2() == t_mat::size1());
	}


	t_mat matRet(mat1.size1(), mat2.size2());
	const std::size_t innersize = std::min(mat1.size2(), mat2.size1());

	for(std::size_t row=0; row<matRet.size1(); ++row)
	{
		for(std::size_t col=0; col<matRet.size2(); ++col)
		{
			matRet(row, col) = 0;
			for(std::size_t i=0; i<innersize; ++i)
				matRet(row, col) += mat1(row, i) * mat2(i, col);
		}
	}

	return matRet;
}


/**
 * element-wise division
 */
template<class t_mat>
t_mat div_perelem(const t_mat& mat1, const t_mat& mat2, bool assert_sizes=true)
requires tl2::is_basic_mat<t_mat>
{
	if(assert_sizes)
	{
		if constexpr(tl2::is_dyn_mat<t_mat>)
			assert(mat1.size1() == mat2.size1() && mat1.size2() == mat2.size2());
	}


	t_mat matRet = zero<t_mat>(mat1.size1(), mat1.size2());

	for(std::size_t row=0; row<matRet.size1(); ++row)
		for(std::size_t col=0; col<matRet.size2(); ++col)
			matRet(row, col) = mat1(row, col) / mat2(row,col);

	return matRet;
}


/**
 * 2-norm
 */
template<class t_vec, class t_real = typename t_vec::value_type>
t_real norm(const t_vec& vec)
requires is_basic_vec<t_vec>
{
	t_real d = static_cast<t_real>(inner<t_vec>(vec, vec));
	return std::sqrt(d);
}


/**
 * n-norm
 */
template<class t_vec, class t_real = typename t_vec::value_type>
typename t_vec::value_type norm(const t_vec& vec, t_real n)
requires is_basic_vec<t_vec>
{
	t_real d = t_real{0};
	for(std::size_t i=0; i<vec.size(); ++i)
		d += std::pow(std::abs(vec[i]), n);
	return std::pow(d, t_real(1)/n);
}


/**
 * outer product
 */
template<class t_mat, class t_vec>
t_mat outer(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	const std::size_t N1 = vec1.size();
	const std::size_t N2 = vec2.size();

	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	for(std::size_t n1=0; n1<N1; ++n1)
	{
		for(std::size_t n2=0; n2<N2; ++n2)
		{
			if constexpr(is_complex<typename t_vec::value_type>)
				mat(n1, n2) = std::conj(vec1[n1]) * vec2[n2];
			else
				mat(n1, n2) = vec1[n1]*vec2[n2];
		}
	}

	return mat;
}


/**
 * outer product (without conjugation of complex vector)
 */
template<class t_mat, class t_vec>
t_mat outer_noconj(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	const std::size_t N1 = vec1.size();
	const std::size_t N2 = vec2.size();

	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	for(std::size_t n1=0; n1<N1; ++n1)
		for(std::size_t n2=0; n2<N2; ++n2)
			mat(n1, n2) = vec1[n1]*vec2[n2];

	return mat;
}



// ----------------------------------------------------------------------------
// operations with metric
// ----------------------------------------------------------------------------

/**
 * covariant metric tensor: g_{i,j} = e_i * e_j
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec, template<class...> class t_cont=std::initializer_list>
t_mat metric(const t_cont<t_vec>& basis_co)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	const std::size_t N = basis_co.size();

	t_mat g_co;
	if constexpr(is_dyn_mat<t_mat>)
		g_co = t_mat(N, N);

	auto iter_i = basis_co.begin();
	for(std::size_t i=0; i<N; ++i)
	{
		auto iter_j = basis_co.begin();
		for(std::size_t j=0; j<N; ++j)
		{
			g_co(i,j) = inner<t_vec>(*iter_i, *iter_j);
			std::advance(iter_j, 1);
		}
		std::advance(iter_i, 1);
	}

	return g_co;
}


/**
 * covariant metric tensor: g_{i,j} = e_i * e_j
 * @see (Arens 2015), p. 815
 */
template<class t_mat>
t_mat metric(const t_mat& basis_co)
requires is_basic_mat<t_mat>
{
	t_mat basis_trans = trans(basis_co);
	return basis_trans * basis_co;
}


/**
 * get levi-civita symbol in fractional coordinates
 * @see (Arens 2015), p. 815
 */
template<class t_mat>
typename t_mat::value_type levi(const t_mat& basis_co,
	const std::initializer_list<std::size_t>& indices)
requires is_basic_mat<t_mat>
{
	using size_t = decltype(basis_co.size1());
	using t_vec = vec<typename t_mat::value_type>;
	t_mat mat = create<t_mat>(basis_co.size1(), basis_co.size2());

	auto iter = indices.begin();
	size_t maxcols = std::min(indices.size(), basis_co.size2());
	for(size_t i=0; i<maxcols; ++i)
	{
		set_col<t_mat, t_vec>(
			mat, col<t_mat, t_vec>(basis_co, *iter), i);
		std::advance(iter, 1);
	}

	return det<t_mat>(mat);
}


/**
 * cross product in fractional coordinates: c^l = eps_ijk g^li a^j b^k
 * @see (Arens 2015), p. 815
 */
template<class t_mat, class t_vec>
t_vec cross(const t_mat& B, const t_vec& a, const t_vec& b)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	using size_t = decltype(B.size1());
	using t_real = typename t_mat::value_type;

	// metric and its inverse
	auto G = metric<t_mat>(B);
	auto [G_inv, ok] = inv(G);

	// maximum indices
	const size_t i_max = G_inv.size1();
	const size_t j_max = a.size();
	const size_t k_max = b.size();
	const size_t l_max = G_inv.size2();

	// cross product result vector
	t_vec c = zero<t_vec>(l_max);

	for(size_t i=0; i<i_max; ++i)
	{
		for(size_t j=0; j<j_max; ++j)
		{
			for(size_t k=0; k<k_max; ++k)
			{
				t_real eps = levi<t_mat>(B, {i,j,k});
				for(size_t l=0; l<l_max; ++l)
					c[l] += eps * G_inv(l,i) * a[j] * b[k];
			}
		}
	}

	return c;
}


/**
 * lower index using metric
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec>
t_vec lower_index(const t_mat& metric_co, const t_vec& vec_contra)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	const std::size_t N = vec_contra.size();
	t_vec vec_co = zero<t_vec>(N);

	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<N; ++j)
			vec_co[i] += metric_co(i,j) * vec_contra[j];

	return vec_co;
}


/**
 * raise index using metric
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec>
t_vec raise_index(const t_mat& metric_contra, const t_vec& vec_co)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	const std::size_t N = vec_co.size();
	t_vec vec_contra = zero<t_vec>(N);

	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<N; ++j)
			vec_contra[i] += metric_contra(i,j) * vec_co[j];

	return vec_contra;
}


/**
 * inner product using metric
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec>
typename t_vec::value_type inner(const t_mat& metric_co,
	const t_vec& vec1_contra, const t_vec& vec2_contra)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	t_vec vec2_co = lower_index<t_mat, t_vec>(metric_co, vec2_contra);
	return inner<t_vec>(vec1_contra, vec2_co);
}


/**
 * 2-norm using metric
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec>
typename t_vec::value_type norm(const t_mat& metric_co, const t_vec& vec_contra)
requires is_basic_vec<t_vec>
{
	return std::sqrt(inner<t_mat, t_vec>(metric_co, vec_contra, vec_contra));
}


/**
 * angle between vectors under a given metric
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec>
typename t_vec::value_type angle(const t_mat& metric_co,
	const t_vec& vec1_contra, const t_vec& vec2_contra)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	using t_real = typename t_mat::value_type;

	t_real len1 = norm<t_mat, t_vec>(metric_co, vec1_contra);
	t_real len2 = norm<t_mat, t_vec>(metric_co, vec2_contra);

	t_real c = inner<t_mat, t_vec>(metric_co, vec1_contra, vec2_contra);
	c /= len1 * len2;
	c = clamp<t_real>(c, -1, 1);

	return std::acos(c);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// projection operators
// ----------------------------------------------------------------------------

/**
 * matrix to project onto vector: P = |v><v|
 * from: |x'> = <v|x> * |v> = |v><v|x> = |v><v| * |x>
 * @see (Arens 2015), p. 814 for the projection tensor
 */
template<class t_mat, class t_vec>
t_mat projector(const t_vec& vec, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	if(is_normalised)
	{
		return outer<t_mat, t_vec>(vec, vec);
	}
	else
	{
		const auto len = norm<t_vec>(vec);
		t_vec _vec = vec / len;
		return outer<t_mat, t_vec>(_vec, _vec);
	}
}


/**
 * project vec1 onto vec2
 *
 * proj_op = |vec2><vec2¦/ len(vec2)^2,  len(vec2) = sqrt(<vec2|vec2>)
 * proj = proj_op * vec1 = |vec2> * <vec2|vec1> / <vec2|vec2>
 *
 * @see (Arens 2015), p. 814 for the projection tensor
 */
template<class t_vec>
t_vec project(const t_vec& vec, const t_vec& vecProj, bool is_normalised = true)
requires is_vec<t_vec>
{
	if(is_normalised)
	{
		return inner<t_vec>(vecProj, vec) * vecProj;
	}
	else
	{
		const auto len = norm<t_vec>(vecProj);
		const t_vec _vecProj = vecProj / len;
		return inner<t_vec>(_vecProj, vec) * _vecProj;
	}
}


/**
 * project vector vec onto another vector vecProj
 * don't multiply with direction vector
 *
 * @see (Arens 2015), p. 814 for the projection tensor
 */
template<class t_vec, class t_real = typename t_vec::value_type>
t_real project_scalar(const t_vec& vec, const t_vec& vecProj, bool is_normalised = true)
requires is_vec<t_vec>
{
	if(is_normalised)
	{
		return inner<t_vec>(vecProj, vec);
	}
	else
	{
		const auto len = norm<t_vec>(vecProj);
		const t_vec _vecProj = vecProj / len;
		return inner<t_vec>(_vecProj, vec);
	}
}


/**
 * project vector vec onto the line lineOrigin + lam*lineDir
 * (shifts line to go through origin, calculate projection and shift back)
 * @returns [closest point, distance, projection parameter]
 *
 * @see https://de.wikipedia.org/wiki/Lot_(Mathematik)
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::tuple<t_vec, t_real, t_real> project_line(const t_vec& vec,
	const t_vec& lineOrigin, const t_vec& _lineDir, bool is_normalised = true)
requires is_vec<t_vec>
{
	const t_real lenDir = is_normalised ? 1 : norm<t_vec>(_lineDir);
	const t_vec lineDir = _lineDir / lenDir;
	const t_vec ptShifted = vec - lineOrigin;

	const t_real paramProj = project_scalar<t_vec, t_real>(ptShifted, lineDir, true);
	const t_vec ptProj = paramProj * lineDir;

	const t_vec ptNearest = lineOrigin + ptProj;
	const t_real dist = norm<t_vec>(vec - ptNearest);
	return std::make_tuple(ptNearest, dist, paramProj);
}


/**
 * distance between point and line
 * @see (Arens 2015), p. 711
 */
template<class t_vec, class t_real = typename t_vec::value_type>
t_real dist_pt_line(const t_vec& pt,
	const t_vec& linePt1, const t_vec& linePt2,
	bool bLineIsInfinite = true)
requires is_vec<t_vec>
{
	const std::size_t dim = linePt1.size();

	const t_vec lineDir = linePt2 - linePt1;
	const auto [nearestPt, dist, paramProj] =
		project_line<t_vec>(pt, linePt1, lineDir, false);


	// get point component with max. difference
	t_real diff = -1.;
	std::size_t compidx = 0;
	for(std::size_t i=0; i<dim; ++i)
	{
		t_real newdiff = std::abs(linePt2[i] - linePt1[i]);
		if(newdiff > diff)
		{
			diff = newdiff;
			compidx = i;
		}
	}


	t_real t = (nearestPt[compidx]-linePt1[compidx]) / (linePt2[compidx]-linePt1[compidx]);
	if(bLineIsInfinite || (t>=t_real{0} && t<=t_real{1}))
	{
		// projection is on line -> use distance between point and projection
		return dist;
	}
	else
	{
		// projection is not on line -> use distance between point and closest line end point
		if(std::abs(t-t_real{0}) < std::abs(t-t_real{1}))
			return norm<t_vec>(linePt1 - pt);
		else
			return norm<t_vec>(linePt2 - pt);
	}
}


/**
 * matrix to project onto orthogonal complement (plane perpendicular to vector): P = 1-|v><v|
 * from completeness relation: 1 = sum_i |v_i><v_i| = |x><x| + |y><y| + |z><z|
 *
 * @see (Arens 2015), p. 814 for the projection tensor
 */
template<class t_mat, class t_vec>
t_mat ortho_projector(const t_vec& vec, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	const std::size_t iSize = vec.size();
	return unit<t_mat>(iSize) -
		projector<t_mat, t_vec>(vec, is_normalised);
}


/**
 * matrix to mirror on plane perpendicular to vector: P = 1 - 2*|v><v|
 * subtracts twice its projection onto the plane normal from the vector
 *
 * @see (Arens 2015), p. 710
 */
template<class t_mat, class t_vec>
t_mat ortho_mirror_op(const t_vec& vec, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using T = typename t_vec::value_type;
	const std::size_t iSize = vec.size();

	return unit<t_mat>(iSize) -
		T(2)*projector<t_mat, t_vec>(vec, is_normalised);
}


/**
 * matrix to mirror [a, b, c, ...] into, e.g.,  [a, b', 0, 0]
 * @see (Scarpino 2011), p. 268
 */
template<class t_mat, class t_vec>
t_mat ortho_mirror_zero_op(const t_vec& vec, std::size_t row)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using T = typename t_vec::value_type;
	const std::size_t N = vec.size();

	t_vec vecSub = zero<t_vec>(N);
	for(std::size_t i=0; i<row; ++i)
		vecSub[i] = vec[i];

	// norm of rest vector
	T n = T(0);
	for(std::size_t i=row; i<N; ++i)
		n += vec[i]*vec[i];
	vecSub[row] = std::sqrt(n);

	const t_vec vecOp = vec - vecSub;

	// nothing to do -> return unit matrix
	if(equals_0<t_vec>(vecOp))
		return unit<t_mat>(vecOp.size(), vecOp.size());

	return ortho_mirror_op<t_mat, t_vec>(vecOp, false);
}


/**
 * project vector vec onto plane through the origin and perpendicular to vector vecNorm
 * (e.g. used to calculate magnetic interaction vector M_perp)
 */
template<class t_vec>
t_vec ortho_project(const t_vec& vec, const t_vec& vecNorm, bool is_normalised = true)
requires is_vec<t_vec>
{
	return vec - project<t_vec>(vec, vecNorm, is_normalised);
}


/**
 * project vector vec onto plane perpendicular to vector vecNorm with distance d
 * vecNorm has to be normalised and plane in Hessian form: x*vecNorm = d
 */
template<class t_vec>
t_vec ortho_project_plane(const t_vec& vec,
	const t_vec& vecNorm, typename t_vec::value_type d)
requires is_vec<t_vec>
{
	// project onto plane through origin
	t_vec vecProj0 = ortho_project<t_vec>(vec, vecNorm, 1);
	// add distance of plane to origin
	return vecProj0 + d*vecNorm;
}


/**
 * mirror a vector on a plane perpendicular to vector vecNorm with distance d
 * vecNorm has to be normalised and plane in Hessian form: x*vecNorm = d
 * @see (Arens 2015), p. 710
 */
template<class t_vec>
t_vec ortho_mirror_plane(const t_vec& vec,
	const t_vec& vecNorm, typename t_vec::value_type d)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec vecProj = ortho_project_plane<t_vec>(vec, vecNorm, d);
	return vec - T(2)*(vec - vecProj);
}


/**
 * find orthonormal substitute basis for vector space (Gram-Schmidt algo)
 * remove orthogonal projections to all other base vectors: |i'> = (1 - sum_{j<i} |j><j|) |i>
 *
 * @see https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
 * @see (Arens 2015), p. 744
 */
template<class t_vec,
	template<class...> class t_cont_in = std::initializer_list,
	template<class...> class t_cont_out = std::vector>
t_cont_out<t_vec> orthonorm_sys(const t_cont_in<t_vec>& sys)
requires is_vec<t_vec>
{
	t_cont_out<t_vec> newsys;
	newsys.reserve(sys.size());

	for(const t_vec& vecSys : sys)
	{
		t_vec vecOrthoProj = vecSys;

		// subtract projections to other basis vectors
		for(const t_vec& vecNewSys : newsys)
			vecOrthoProj -= project<t_vec>(vecSys, vecNewSys, true);

		// normalise
		vecOrthoProj /= norm<t_vec>(vecOrthoProj);
		newsys.emplace_back(std::move(vecOrthoProj));
	}

	return newsys;
}
// ----------------------------------------------------------------------------


/**
 * submatrix removing a column/row from a matrix stored in a vector container
 */
template<class t_vec, class t_matvec = t_vec>
t_vec flat_submat(const t_matvec& mat,
	std::size_t iNumRows, std::size_t iNumCols,
	std::size_t iRemRow, std::size_t iRemCol)
requires is_basic_vec<t_vec>
{
	t_vec vec;
	vec.reserve(mat.size());

	for(std::size_t iRow=0; iRow<iNumRows; ++iRow)
	{
		if(iRow == iRemRow)
			continue;

		for(std::size_t iCol=0; iCol<iNumCols; ++iCol)
		{
			if(iCol == iRemCol)
				continue;
			vec.push_back(mat[iRow*iNumCols + iCol]);
		}
	}

	return vec;
}


/**
 * submatrix removing a column/row from a matrix
 */
template<class t_mat>
t_mat submat(const t_mat& mat, decltype(mat.size1()) iRemRow, decltype(mat.size2()) iRemCol)
requires is_basic_mat<t_mat> && is_dyn_mat<t_mat>
{
	using size_t = decltype(mat.size1());
	t_mat matRet = create<t_mat>(mat.size1()-1, mat.size2()-1);

	size_t iResRow = 0;
	for(size_t iRow=0; iRow<mat.size1(); ++iRow)
	{
		if(iRow == iRemRow)
			continue;

		size_t iResCol = 0;
		for(size_t iCol=0; iCol<mat.size2(); ++iCol)
		{
			if(iCol == iRemCol)
				continue;

			matRet(iResRow, iResCol) = mat(iRow, iCol);
			++iResCol;
		}

		++iResRow;
	}

	return matRet;
}


/**
 * submatrix of a given size starting at given indices
 */
template<class t_mat>
t_mat submat(const t_mat& mat,
	decltype(mat.size1()) row_start, decltype(mat.size2()) col_start,
	decltype(mat.size1()) num_rows, decltype(mat.size2()) num_cols)
requires is_basic_mat<t_mat> && is_dyn_mat<t_mat>
{
	using size_t = decltype(mat.size1());
	t_mat matRet = create<t_mat>(num_rows, num_cols);

	for(size_t row=0; row<num_rows; ++row)
		for(size_t col=0; col<num_cols; ++col)
			matRet(row, col) = mat(row+row_start, col+col_start);

	return matRet;
}


/**
 * set submatrix at given starting indices
 */
template<class t_mat>
void set_submat(t_mat& mat, const t_mat& submat,
	decltype(mat.size1()) row_start, decltype(mat.size2()) col_start,
	decltype(mat.size1()) num_rows=9999, decltype(mat.size2()) num_cols=9999)
requires is_basic_mat<t_mat> && is_dyn_mat<t_mat>
{
	using size_t = decltype(mat.size1());

	num_rows = std::min(num_rows, submat.size1());
	num_cols = std::min(num_cols, submat.size2());

	for(size_t row=0; row<num_rows; ++row)
		for(size_t col=0; col<num_cols; ++col)
			mat(row+row_start, col+col_start) = submat(row, col);
}


/**
 * add submatrix at given starting indices
 */
template<class t_mat>
void add_submat(t_mat& mat, const t_mat& submat,
	decltype(mat.size1()) row_start, decltype(mat.size2()) col_start)
requires is_basic_mat<t_mat> && is_dyn_mat<t_mat>
{
	using size_t = decltype(mat.size1());

	for(size_t row=0; row<submat.size1(); ++row)
	{
		for(size_t col=0; col<submat.size2(); ++col)
			mat(row+row_start, col+col_start) += submat(row, col);
	}
}


/**
 * determinant from a square matrix stored in a vector container
 * @see (Merziger 2006), p. 185
 */
template<class t_vec, class t_matvec = t_vec>
typename t_vec::value_type flat_det(const t_matvec& mat, std::size_t iN)
requires is_basic_vec<t_vec>
{
	using T = typename t_vec::value_type;

	// special cases
	if(iN == 0)
		return 0;
	else if(iN == 1)
		return mat[0];
	else if(iN == 2)
		return mat[0]*mat[3] - mat[1]*mat[2];

	// recursively expand determiant along a row
	T fullDet = T(0);
	std::size_t iRow = 0;

	// get row with maximum number of zeros
	std::size_t iMaxNumZeros = 0;
	for(std::size_t iCurRow=0; iCurRow<iN; ++iCurRow)
	{
		std::size_t iNumZeros = 0;
		for(std::size_t iCurCol=0; iCurCol<iN; ++iCurCol)
		{
			if(equals<T>(mat[iCurRow*iN + iCurCol], T(0)))
				++iNumZeros;
		}

		if(iNumZeros > iMaxNumZeros)
		{
			iRow = iCurRow;
			iMaxNumZeros = iNumZeros;
		}
	}

	for(std::size_t iCol=0; iCol<iN; ++iCol)
	{
		const T elem = mat[iRow*iN + iCol];
		if(equals<T>(elem, 0))
			continue;

		const T sgn = ((iRow+iCol) % 2) == 0 ? T(1) : T(-1);
		const t_vec subMat = flat_submat<t_vec, t_matvec>(
			mat, iN, iN, iRow, iCol);
		const T subDet = flat_det<t_vec>(subMat, iN-1) * sgn;

		fullDet += elem * subDet;
	}

	return fullDet;
}


/**
 * trace
 */
template<class t_mat>
typename t_mat::value_type trace(const t_mat& mat)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;
	T _tr = T(0);

	std::size_t N = std::min(mat.size1(), mat.size2());
	for(std::size_t i=0; i<N; ++i)
		_tr += mat(i,i);

	return _tr;
}


/**
 * gets reciprocal basis vectors |b_i> from real basis vectors |a_i> (and vice versa)
 * c: multiplicative constant (c=2*pi for physical lattices, c=1 for mathematics)
 *
 * Def.: <b_i | a_j> = c * delta(i,j)  =>
 *
 * e.g. 2d case:
 *                   ( a_1x  a_2x )
 *                   ( a_1y  a_2y )
 *
 * ( b_1x  b_1y )    (    1     0 )
 * ( b_2x  b_2y )    (    0     1 )
 *
 * B^t * A = I
 * A = B^(-t)
 */
template<class t_mat, class t_vec,
	template<class...> class t_cont_in = std::initializer_list,
	template<class...> class t_cont_out = std::vector>
t_cont_out<t_vec> recip(const t_cont_in<t_vec>& lstReal, typename t_vec::value_type c=1)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	const t_mat basis = create<t_mat, t_vec, t_cont_in>(lstReal);
	auto [basis_inv, bOk] = inv<t_mat>(basis);
	basis_inv *= c;

	t_cont_out<t_vec> lstRecip;
	lstRecip.reserve(basis_inv.size1());

	for(std::size_t currow=0; currow<basis_inv.size1(); ++currow)
	{
		const t_vec rowvec = row<t_mat, t_vec>(basis_inv, currow);
		lstRecip.emplace_back(std::move(rowvec));
	}

	return lstRecip;
}


template<class t_vec>
t_vec cross(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec>;


/**
 * general n-dim cross product using determinant definition
 * @see https://en.wikipedia.org/wiki/Cross_product
 */
template<class t_vec, template<class...> class t_cont = std::initializer_list>
t_vec cross(const t_cont<t_vec>& vecs)
requires is_basic_vec<t_vec>
{
	using t_size = decltype(t_vec{}.size());
	using T = typename t_vec::value_type;

	// N also has to be equal to the vector size!
	const t_size N = vecs.size()+1;
	t_vec vec = zero<t_vec>(N);

	// 3-dim case
	if(N == 3 && vecs.begin()->size() == 3)
	{
		return cross<t_vec>(*vecs.begin(), *std::next(vecs.begin(), 1));
	}

	// general case
	else
	{
		for(t_size iComp=0; iComp<N; ++iComp)
		{
			std::vector<T> mat = zero<std::vector<T>>(N*N);
			mat[0*N + iComp] = T(1);

			t_size row_idx = 0;
			for(const t_vec& vec : vecs)
			{
				for(t_size col_idx=0; col_idx<N; ++col_idx)
					mat[(row_idx+1)*N + col_idx] = vec[col_idx];
				++row_idx;
			}

			vec[iComp] = flat_det<decltype(mat)>(mat, N);
		}
	}

	return vec;
}


/**
 * intersection of plane <x|n> = d and line |org> + lam*|dir>
 * @returns [position of intersection, 0: no intersection, 1: intersection, 2: line on plane, line parameter lambda]
 *
 * insert |x> = |org> + lam*|dir> in plane equation:
 * <org|n> + lam*<dir|n> = d
 * lam = (d - <org|n>) / <dir|n>
 *
 * @see http://mathworld.wolfram.com/Line-PlaneIntersection.html
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::tuple<t_vec, int, t_real>
intersect_line_plane(
	const t_vec& lineOrg, const t_vec& lineDir,
	const t_vec& planeNorm, t_real plane_d,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_vec<t_vec>
{
	// are line and plane parallel?
	const t_real dir_n = inner<t_vec>(lineDir, planeNorm);
	if(equals<t_real>(dir_n, 0, eps))
	{
		const t_real org_n = inner<t_vec>(lineOrg, planeNorm);
		// line on plane?
		if(equals<t_real>(org_n, plane_d, eps))
			return std::make_tuple(t_vec(), 2, t_real(0));
		// no intersection
		return std::make_tuple(t_vec(), 0, t_real(0));
	}

	const t_real org_n = inner<t_vec>(lineOrg, planeNorm);
	const t_real lam = (plane_d - org_n) / dir_n;

	const t_vec vecInters = lineOrg + lam*lineDir;
	return std::make_tuple(vecInters, 1, lam);
}


/**
 * intersection of plane <x|n> = d and a polygon
 * @return vertices of plane - polygon edge intersections
 */
template<class t_vec, class t_real = typename t_vec::value_type,
	template<class...> class t_cont>
t_cont<t_vec> intersect_plane_poly(
	const t_vec& planeNorm, t_real plane_d,
	const t_cont<t_vec>& polyVerts,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_vec<t_vec>
{
	t_cont<t_vec> edgeInters;

	// intersect with each polygon edge
	for(std::size_t i = 0; i < polyVerts.size(); ++i)
	{
		std::size_t j = (i+1) % polyVerts.size();

		t_vec lineOrg = polyVerts[i];
		t_vec lineDir = polyVerts[j] - lineOrg;

		auto [pos, ty, lam] = intersect_line_plane(
			lineOrg, lineDir, planeNorm, plane_d, eps);

		if(ty == 1 && lam >= 0. && lam < 1.)
		{
			edgeInters.emplace_back(std::move(pos));
		}
		if(ty == 2)
		{
			edgeInters.push_back(polyVerts[i]);
			edgeInters.push_back(polyVerts[j]);
		}
	}

	return edgeInters;
}


/**
 * intersection of a sphere and a line |org> + lam*|dir>
 * @returns vector of intersections
 * insert |x> = |org> + lam*|dir> in sphere equation <x-mid | x-mid> = r^2
 *
 * @see https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection for the solution.
 */
template<class t_vec, template<class...> class t_cont = std::vector>
requires is_vec<t_vec>
t_cont<t_vec> intersect_line_sphere(
	const t_vec& lineOrg, const t_vec& _lineDir,
	const t_vec& sphereOrg, typename t_vec::value_type sphereRad,
	bool linedir_normalised = false, bool only_segment = false,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
{
	using T = typename t_vec::value_type;

	t_vec lineDir = _lineDir;
	T lenDir = linedir_normalised ? T(1) : norm<t_vec>(lineDir);

	if(!linedir_normalised)
		lineDir /= lenDir;

	auto vecDiff = sphereOrg - lineOrg;
	auto proj = project_scalar<t_vec>(vecDiff, lineDir, true);
	auto rt = proj*proj + sphereRad*sphereRad - inner<t_vec>(vecDiff, vecDiff);

	// no intersection
	if(rt < T(0))
		return t_cont<t_vec>{};

	// one intersection
	if(equals(rt, T(0), eps))
	{
		T lam = proj/lenDir;
		if(!only_segment || (only_segment && lam >= T(0) && lam < T(1)))
			return t_cont<t_vec>{{ lineOrg + proj*lineDir }};
		return t_cont<t_vec>{};
	}

	// two intersections
	auto val = std::sqrt(rt);
	t_cont<t_vec> inters;
	inters.reserve(2);

	T lam1 = (proj + val)/lenDir;
	T lam2 = (proj - val)/lenDir;

	if(!only_segment || (only_segment && lam1 >= T(0) && lam1 < T(1)))
		inters.emplace_back(lineOrg + (proj + val)*lineDir);
	if(!only_segment || (only_segment && lam2 >= T(0) && lam2 < T(1)))
		inters.emplace_back(lineOrg + (proj - val)*lineDir);

	// sort intersections by x
	std::sort(inters.begin(), inters.end(),
		[](const t_vec& vec1, const t_vec& vec2) -> bool
		{ return vec1[0] < vec2[0]; });

	return inters;
}


/**
 * average vector or matrix
 */
template<class ty, template<class...> class t_cont = std::vector>
ty avg(const t_cont<ty>& vecs)
requires is_vec<ty> || is_mat<ty>
{
	if(vecs.size() == 0)
		return ty();

	typename ty::value_type num = 1;
	ty vec = *vecs.begin();

	auto iter = vecs.begin();
	std::advance(iter, 1);

	for(; iter!=vecs.end(); std::advance(iter, 1))
	{
		vec += *iter;
		++num;
	}
	vec /= num;

	return vec;
}


/**
 * intersection of a polygon and a line
 * @returns [position of intersection, intersects?, line parameter lambda]
 */
template<class t_vec, template<class ...> class t_cont = std::vector>
std::tuple<t_vec, bool, typename t_vec::value_type>
intersect_line_poly(
	const t_vec& lineOrg, const t_vec& lineDir,
	const t_cont<t_vec>& poly)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	// middle point
	const t_vec mid = avg<t_vec, t_cont>(poly);

	// calculate polygon plane
	const t_vec vec0 = poly[0] - mid;
	const t_vec vec1 = poly[1] - mid;
	t_vec planeNorm = cross<t_vec>({vec0, vec1});
	planeNorm /= norm<t_vec>(planeNorm);
	const T planeD = inner<t_vec>(poly[0], planeNorm);

	// intersection with plane
	auto [vec, intersects, lam] =
		intersect_line_plane<t_vec>(lineOrg, lineDir, planeNorm, planeD);
	if(intersects != 1)
		return std::make_tuple(t_vec(), false, T(0));

	// is intersection point contained in polygon?
	const t_vec* vecFirst = &(*poly.rbegin());
	for(auto iter=poly.begin(); iter!=poly.end(); std::advance(iter, 1))
	{
		const t_vec* vecSecond = &(*iter);
		const t_vec edge = *vecSecond - *vecFirst;

		// plane through edge
		t_vec edgeNorm = cross<t_vec>({edge, planeNorm});
		edgeNorm /= norm<t_vec>(edgeNorm);
		const T edgePlaneD = inner<t_vec>(*vecFirst, edgeNorm);

		// side of intersection
		const T ptEdgeD = inner<t_vec>(vec, edgeNorm);

		// outside polygon?
		if(ptEdgeD > edgePlaneD)
			return std::make_tuple(t_vec(), false, T(0));

		vecFirst = vecSecond;
	}

	// intersects with polygon
	return std::make_tuple(vec, true, lam);
}


/**
 * intersection of a polygon (transformed with a matrix) and a line
 * @returns [position of intersection, intersects?, line parameter lambda]
 */
template<class t_vec, class t_mat, template<class ...> class t_cont = std::vector>
std::tuple<t_vec, bool, typename t_vec::value_type>
intersect_line_poly(
	const t_vec& lineOrg, const t_vec& lineDir,
	const t_cont<t_vec>& _poly, const t_mat& mat)
requires is_vec<t_vec> && is_mat<t_mat>
{
	auto poly = _poly;

	// transform each vertex of the polygon
	// TODO: check for homogeneous coordinates!
	for(t_vec& vec : poly)
		vec = mat * vec;

	return intersect_line_poly<t_vec, t_cont>(lineOrg, lineDir, poly);
}


/**
 * intersection or closest points of lines |org1> + lam1*|dir1> and |org2> + lam2*|dir2>
 * @returns [nearest position 1, nearest position 2, dist, valid?, line parameter 1, line parameter 2]
 *
 * |org1> + lam1*|dir1>  =  |org2> + lam2*|dir2>
 * |org1> - |org2>  =  lam2*|dir2> - lam1*|dir1>
 * |org1> - |org2>  =  (dir2 | -dir1) * |lam2 lam1>
 * (dir2 | -dir1)^T * (|org1> - |org2>)  =  (dir2 | -dir1)^T * (dir2 | -dir1) * |lam2 lam1>
 * |lam2 lam1> = ((dir2 | -dir1)^T * (dir2 | -dir1))^(-1) * (dir2 | -dir1)^T * (|org1> - |org2>)
 *
 * @see https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::tuple<t_vec, t_vec, bool, t_real, t_real, t_real>
intersect_line_line(
	const t_vec& line1Org, const t_vec& line1Dir,
	const t_vec& line2Org, const t_vec& line2Dir,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_vec<t_vec>
{
	const t_vec orgdiff = line1Org - line2Org;

	// direction matrix (symmetric)
	const t_real d11 = inner<t_vec>(line2Dir, line2Dir);
	const t_real d12 = -inner<t_vec>(line2Dir, line1Dir);
	const t_real d22 = inner<t_vec>(line1Dir, line1Dir);

	const t_real d_det = d11*d22 - d12*d12;

	// check if matrix is invertible
	if(equals<t_real>(d_det, 0, eps))
		return std::make_tuple(t_vec(), t_vec(), false, 0, 0, 0);

	// inverse (symmetric)
	const t_real d11_i = d22 / d_det;
	const t_real d12_i = -d12 / d_det;
	const t_real d22_i = d11 / d_det;

	const t_vec v1 = d11_i*line2Dir - d12_i*line1Dir;
	const t_vec v2 = d12_i*line2Dir - d22_i*line1Dir;

	const t_real lam2 = inner<t_vec>(v1, orgdiff);
	const t_real lam1 = inner<t_vec>(v2, orgdiff);

	const t_vec pos1 = line1Org + lam1*line1Dir;
	const t_vec pos2 = line2Org + lam2*line2Dir;
	const t_real dist = norm<t_vec>(pos2 - pos1);

	return std::make_tuple(pos1, pos2, true, dist, lam1, lam2);
}


/**
 * intersection of planes <x|n1> = d1 and <x|n2> = d2
 * @returns line [org, dir, 0: no intersection, 1: intersection, 2: planes coincide]
 *
 * @see http://mathworld.wolfram.com/Plane-PlaneIntersection.html
 */
template<class t_vec>
std::tuple<t_vec, t_vec, int>
	intersect_plane_plane(
	const t_vec& plane1Norm, typename t_vec::value_type plane1_d,
	const t_vec& plane2Norm, typename t_vec::value_type plane2_d)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec lineDir = cross<t_vec>({plane1Norm, plane2Norm});
	const T lenCross = norm<t_vec>(lineDir);

	// planes parallel or coinciding
	if(equals<T>(lenCross, 0))
	{
		const bool bCoincide = equals<T>(plane1_d, plane2_d);
		return std::make_tuple(t_vec(), t_vec(), bCoincide ? 2 : 0);
	}

	lineDir /= lenCross;

	t_vec lineOrg = - cross<t_vec>({plane1Norm, lineDir}) * plane2_d
		+ cross<t_vec>({plane2Norm, lineDir}) * plane1_d;
	lineOrg /= lenCross;

	return std::make_tuple(lineOrg, lineDir, 1);
}


/**
 * uv coordinates of a point inside a polygon defined by three vertices
 */
template<class t_vec>
t_vec poly_uv_ortho(const t_vec& vert1, const t_vec& vert2, const t_vec& vert3,
	const t_vec& uv1, const t_vec& uv2, const t_vec& uv3,
	const t_vec& _pt)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec vec12 = vert2 - vert1;
	t_vec vec13 = vert3 - vert1;

	t_vec uv12 = uv2 - uv1;
	t_vec uv13 = uv3 - uv1;


	// ----------------------------------------------------
	// orthonormalisation
	const T len12 = norm<t_vec>(vec12);
	const T len13 = norm<t_vec>(vec13);
	const T lenuv12 = norm<t_vec>(uv12);
	const T lenuv13 = norm<t_vec>(uv13);
	auto vecBasis = orthonorm_sys<t_vec, std::initializer_list, std::vector>({vec12, vec13});
	auto uvBasis = orthonorm_sys<t_vec, std::initializer_list, std::vector>({uv12, uv13});
	vec12 = vecBasis[0]*len12; vec13 = vecBasis[1]*len13;
	uv12 = uvBasis[0]*lenuv12; uv13 = uvBasis[1]*lenuv13;
	// ----------------------------------------------------


	const t_vec pt = _pt - vert1;

	// project a point onto a vector and return the fraction along that vector
	auto project_lam = [](const t_vec& vec, const t_vec& vecProj) -> T
	{
		const T len = norm<t_vec>(vecProj);
		const t_vec _vecProj = vecProj / len;
		T lam = inner<t_vec>(_vecProj, vec);
		return lam / len;
	};

	T lam12 = project_lam(pt, vec12);
	T lam13 = project_lam(pt, vec13);

	// uv coordinates at specified point
	const t_vec uv_pt = uv1 + lam12*uv12 + lam13*uv13;
	return uv_pt;
}


/**
 * uv coordinates of a point inside a polygon defined by three vertices
 * (more general version than poly_uv_ortho)
 */
template<class t_mat, class t_vec>
t_vec poly_uv(const t_vec& vert1, const t_vec& vert2, const t_vec& vert3,
	const t_vec& uv1, const t_vec& uv2, const t_vec& uv3,
	const t_vec& _pt)
requires is_mat<t_mat> && is_vec<t_vec>
{
	t_vec vec12 = vert2 - vert1;
	t_vec vec13 = vert3 - vert1;
	t_vec vecnorm = cross<t_vec>({vec12, vec13});

	// basis
	const t_mat basis = create<t_mat, t_vec>({vec12, vec13, vecnorm}, false);

	// reciprocal basis, RECI = REAL^(-T)
	const auto [basisInv, bOk] = inv<t_mat>(basis);
	if(!bOk) return zero<t_vec>(uv1.size());

	t_vec pt = _pt - vert1;		// real pt
	pt = basisInv * pt;		// reciprocal pt

	// uv coordinates at specified point
	t_vec uv12 = uv2 - uv1;
	t_vec uv13 = uv3 - uv1;

	// pt has components in common reciprocal basis
	// assumes that both vector and uv coordinates have the same reciprocal basis
	return uv1 + pt[0]*uv12 + pt[1]*uv13;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// statistical functions
// ----------------------------------------------------------------------------

/**
 * mean value
 */
template<class t_elem, template<class...> class t_cont = std::vector>
t_elem mean(const t_cont<t_elem>& vec)
requires is_basic_vec<t_cont<t_elem>>
{
	if(vec.size()==0) return t_elem{};
	else if(vec.size()==1) return *vec.begin();

	using namespace tl2_ops;

	//t_elem meanvec = std::accumulate(std::next(vec.begin(), 1), vec.end(), *vec.begin());

	t_elem meanvec = *vec.begin();
	auto iter = std::next(vec.begin(), 1);
	for(; iter!=vec.end(); iter=std::next(iter, 1))
		meanvec += *iter;

	meanvec /= vec.size();
	return meanvec;
}


/**
 * mean value with given probability
 */
template<class t_vec_prob, class t_vec>
typename t_vec::value_type mean(const t_vec_prob& vecP, const t_vec& vec)
requires is_basic_vec<t_vec> && is_basic_vec<t_vec_prob>
{
	typedef typename t_vec::value_type T;
	typedef typename t_vec_prob::value_type Tprob;
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
 * @see https://en.wikipedia.org/wiki/Bessel%27s_correction
 */
template<class t_vec>
typename t_vec::value_type std_dev(const t_vec& vec, bool bCorr=1)
requires is_basic_vec<t_vec>
{
	typedef typename t_vec::value_type T;
	if(vec.size()<=1) return T(0);

	T tProb = T(vec.size());
	if(bCorr) tProb -= T(1);

	T tMean = mean(vec);
	T t = T(0);
	for(const T& tval : vec)
		t += (tval-tMean) * (tval-tMean);
	t /= tProb;

	return std::sqrt(t);
}


/**
 * standard deviation with given probability
 * @see https://en.wikipedia.org/wiki/Standard_deviation
 */
template<class t_vec_prob, class t_vec>
typename t_vec::value_type std_dev(const t_vec_prob& vecP, const t_vec& vec)
requires is_basic_vec<t_vec> && is_basic_vec<t_vec_prob>
{
	typedef typename t_vec::value_type T;
	std::size_t iSize = std::min(vecP.size(), vec.size());
	if(iSize<=1) return T(0);

	T tMean = mean<t_vec_prob, t_vec>(vecP, vec);
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
 * calculates minimum and maximum components of a collection of vectors
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_vec, t_vec> minmax(const t_cont<t_vec>& verts)
requires is_vec<t_vec>
{
	using namespace tl2_ops;
	using t_real = typename t_vec::value_type;

	if(!verts.size())
		return std::make_tuple(t_vec{}, t_vec{});

	// set to limit values
	t_vec vecmin = zero<t_vec>(verts.begin()->size());
	t_vec vecmax = zero<t_vec>(verts.begin()->size());

	for(std::size_t i=0; i<vecmin.size(); ++i)
	{
		vecmin[i] = std::numeric_limits<t_real>::max();
		vecmax[i] = std::numeric_limits<t_real>::lowest();
	}

	// iterate components
	for(std::size_t i=0; i<vecmin.size(); ++i)
	{
		// iterate vectors
		for(const t_vec& vec : verts)
		{
			vecmin[i] = std::min(vecmin[i], vec[i]);
			vecmax[i] = std::max(vecmax[i], vec[i]);
		}
	}

	return std::make_tuple(vecmin, vecmax);
}


/**
 * calculates the covariance and the correlation matrices
 * covariance: C_ij = cov(X_i, X_j) = < (X_i - <X_i>) * (X_j - <X_j>) >
 * correlation: K_ij = C_ij / (sigma_i sigma_j)
 *
 * @see http://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm
 * @see (Arfken 2013) pp. 1142-1144
 * @see (Arens 2015), p. 795 and p. 1372
 */
template<class t_mat, class t_vec, class T=typename t_vec::value_type>
std::tuple<t_mat, t_mat>
covariance(const std::vector<t_vec>& vecVals, const std::vector<T>* pProb = 0)
requires is_mat<t_mat> && is_vec<t_vec>
{
	using t_vecvec = typename std::remove_reference<decltype(vecVals)>::type;
	using t_innervec_org = decltype(vecVals[0]);
	using t_innervec = typename std::remove_const<
		typename std::remove_reference<t_innervec_org>::type>::type;

	if(vecVals.size() == 0)
		return std::make_tuple(t_mat(), t_mat());

	// mean vector <X_i>
	t_innervec vecMean;
	if(pProb)
		vecMean = mean<std::vector<T>, t_vecvec>(*pProb, vecVals);
	else
		vecMean = mean<t_vec>(vecVals);

	t_mat matCov = zero<t_mat>(vecVals[0].size(), vecVals[0].size());
	T tSum = T{0};
	const std::size_t N = vecVals.size();

	for(std::size_t i=0; i<N; ++i)
	{
		T tprob = T{1};

		// X_i - <X_i>
		t_innervec vec = vecVals[i] - vecMean;

		// matrix elements, AA^t
		t_mat matOuter = outer<t_mat, t_vec>(vec, vec);

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
	t_innervec vecVar = diag_vec<t_vec, t_mat>(matCov);
	t_innervec vecStdDev(vecVar.size());

	std::transform(vecVar.begin(), vecVar.end(), vecStdDev.begin(),
		[](typename t_innervec::value_type d) -> typename t_innervec::value_type
		{ return std::sqrt(d); });

	t_mat matStdDev = outer<t_mat, t_vec>(vecStdDev, vecStdDev);
	t_mat matCorr = div_perelem<t_mat>(matCov, matStdDev);
	// --------------------------------------------------------------------------------

	return std::make_tuple(matCov, matCorr);
}


/**
 * calculates chi^2 distance of a function model to data points
 * chi^2 = sum( (y_i - f(x_i))^2 / sigma_i^2 )
 *
 * @see (Arfken 2013), p. 1170
 */
template<class T, class t_func, class t_iter_dat=T*>
T chi2(const t_func& func, std::size_t N,
	const t_iter_dat x, const t_iter_dat y, const t_iter_dat dy)
{
	using t_dat = typename std::remove_pointer<t_iter_dat>::type;
	T tchi2 = T{0};

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


/**
 * chi^2 for vector types
 *
 * @see (Merziger 2006), p. 185
 */
template<class t_vec, class t_func>
typename t_vec::value_type chi2(const t_func& func,
	const t_vec& x, const t_vec& y, const t_vec& dy)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;
	return chi2<T, t_func, T*>(func, x.size(), x.data(), y.data(),
		dy.size() ? dy.data() : nullptr);
}


/**
 * chi^2 which doesn't use an x value, but an index instead: y[idx] - func(idx)
 * @see (Arfken 2013), p. 1170
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
 * @see (Arfken 2013), p. 1170
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
 * @see (Arfken 2013), p. 1170
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
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// 3-dim algos
// ----------------------------------------------------------------------------

/**
 * 3-dim cross product
 * @see https://en.wikipedia.org/wiki/Cross_product
 */
template<class t_vec>
t_vec cross(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec>
{
	t_vec vec;

	// only valid for 3-vectors -> use first three components
	if(vec1.size() < 3 || vec2.size() < 3)
		return vec;

	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(3);

	for(int i=0; i<3; ++i)
		vec[i] = vec1[(i+1)%3]*vec2[(i+2)%3] - vec1[(i+2)%3]*vec2[(i+1)%3];

	return vec;
}


/**
 * cross product matrix (3x3)
 * @see https://en.wikipedia.org/wiki/Skew-symmetric_matrix
 */
template<class t_mat, class t_vec>
t_mat skewsymmetric(const t_vec& vec)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(3,3);

	// if static matrix is larger than 3x3 (e.g. for homogeneous coordinates), initialise as identity
	if(mat.size1() > 3 || mat.size2() > 3)
		mat = unit<t_mat>(mat.size1(), mat.size2());

	mat(0,0) = 0;       mat(0,1) = -vec[2]; mat(0,2) = vec[1];
	mat(1,0) = vec[2];  mat(1,1) = 0;       mat(1,2) = -vec[0];
	mat(2,0) = -vec[1]; mat(2,1) = vec[0];  mat(2,2) = 0;

	return mat;
}


/**
 * givens rotation matrix
 * @see https://en.wikipedia.org/wiki/Givens_rotation
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat givens(std::size_t N, std::size_t i, std::size_t j, t_real angle)
requires is_mat<t_mat>
{
	t_mat mat = unit<t_mat>(N, N);

	const t_real s = std::sin(angle);
	const t_real c = std::cos(angle);

	mat(i, j) = -s;
	mat(j, i) = +s;
	mat(i, i) = mat(j, j) = c;

	return mat;
}


/**
 * givens rotation matrix with precalculated sin and cos constants
 * @see https://en.wikipedia.org/wiki/Givens_rotation
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat givens(std::size_t N, std::size_t i, std::size_t j, t_real s, t_real c)
requires is_mat<t_mat>
{
	t_mat mat = unit<t_mat>(N, N);

	mat(i, j) = -s;
	mat(j, i) = +s;
	mat(i, i) = mat(j, j) = c;

	return mat;
}


/**
 * SO(2) rotation matrix
 * @see https://en.wikipedia.org/wiki/Rotation_matrix
 */
template<class t_mat>
t_mat rotation_2d(const typename t_mat::value_type angle)
requires tl2::is_mat<t_mat>
{
	return givens<t_mat>(2, 0, 1, angle);
}


/**
 * SO(3) matrix to rotate around an axis (Rodrigues' formula)
 * @see (Arens 2015), p. 718 and p. 816
 * @see (Merziger 2006), p. 208
 * @see https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
 */
template<class t_mat, class t_vec>
t_mat rotation(const t_vec& axis, const typename t_vec::value_type angle,
	bool is_normalised = 1)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_real = typename t_vec::value_type;

	const t_real c = std::cos(angle);
	const t_real s = std::sin(angle);

	t_real len = 1;
	if(!is_normalised)
		len = norm<t_vec>(axis);

	// ----------------------------------------------------
	// special cases: rotations around [100], [010], [001]
	if(equals(axis, create<t_vec>({ len, 0, 0 })))
		return create<t_mat>({{1,0,0}, {0,c,s}, {0,-s,c}});
	else if(equals(axis, create<t_vec>({ 0, len, 0 })))
		return create<t_mat>({{c,0,-s}, {0,1,0}, {s,0,c}});
	else if(equals(axis, create<t_vec>({ 0, 0, len })))
		return create<t_mat>({{c,s,0}, {-s,c,0}, {0,0,1}});

	// ----------------------------------------------------
	// general case
	// project along rotation axis
	t_mat matProj1 = projector<t_mat, t_vec>(axis, is_normalised);

	// project along axis 2 in plane perpendicular to rotation axis
	t_mat matProj2 = ortho_projector<t_mat, t_vec>(axis, is_normalised) * c;

	// project along axis 3 in plane perpendicular to rotation axis and axis 2
	t_mat matProj3 = skewsymmetric<t_mat, t_vec>(axis/len) * s;

	//std::cout << matProj1(3,3) <<  " " << matProj2(3,3) <<  " " << matProj3(3,3) << std::endl;
	t_mat matProj = matProj1 + matProj2 + matProj3;

	// if matrix is larger than 3x3 (e.g. for homogeneous cooridnates), fill up with identity
	unit<t_mat>(matProj, 3,3, matProj.size1(), matProj.size2());
	return matProj;
}


/**
 * get a vector perpendicular to vec
 */
template<class t_vec, class t_scalar = typename t_vec::value_type>
requires is_vec<t_vec>
t_vec perp(const t_vec& vec, t_scalar eps = std::numeric_limits<t_scalar>::epsilon())
{
	if(vec.size() == 2)
	{
		t_vec vecret = create<t_vec>({vec[1], vec[0]});
		if(std::abs(vecret[0]) > std::abs(vecret[1]))
			vecret[0] = -vecret[0];
	}

	else if(vec.size() == 3 || vec.size() == 4)
	{
		while(1)
		{
			t_vec rand = tl2::rand<t_vec>(3);
			t_vec perp = cross<t_vec>({ vec, rand });
			t_scalar dot = inner<t_vec>(perp, perp);

			if(dot > eps)
			{
				perp /= std::sqrt(dot);
				return perp;
			}
		}
	}

	return t_vec{};
}


/**
 * matrix to rotate a vector into [1, 0, 0, ...]
 * @see (Zhelezov 2017) O. I. Zhelezov, American Journal of Computational and Applied Mathematics 7(2), pp. 51-57 (2017), doi: 10.5923/j.ajcam.20170702.04
 */
template<class t_mat, class t_vec, class t_real = typename t_vec::value_type>
t_mat rotation_x0(const t_vec& _vec, t_real eps = 1e-6)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_size = decltype(_vec.size());
	const t_size N = _vec.size();

	t_vec vec = _vec;
	t_mat mat = unit<t_mat>(N);

	// equations (6) and (7) from (Zhelezov 2017)
	for(t_size n=1; n<N; ++n)
	{
		t_size i = N - n - 1;
		t_size j = i + 1;

		t_real len = vec[i]*vec[i] + vec[j]*vec[j];
		if(!equals_0<t_real>(len, eps))
		{
			len = std::sqrt(len);
			t_real s = -vec[j] / len;
			t_real c = vec[i] / len;

			t_mat rot_ij = givens<t_mat, t_real>(N, i, j, s, c);
			vec = rot_ij * vec;
			mat = rot_ij * mat;
		}
	}

	return mat;
}


/**
 * matrix to rotate vector vec1 into vec2
 * (the vectors do not need to be normalised)
 * @see (Zhelezov 2017) O. I. Zhelezov, American Journal of Computational and Applied Mathematics 7(2), pp. 51-57 (2017), doi: 10.5923/j.ajcam.20170702.04
 */
template<class t_mat, class t_vec, class t_real = typename t_vec::value_type>
t_mat rotation(const t_vec& vec1, const t_vec& vec2,
	const t_vec *perp_vec = nullptr, t_real eps = 1e-6, bool force_3dim = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	assert(vec1.size() == vec2.size());

	using t_size = decltype(vec1.size());
	const t_size dim = vec1.size();

	// 2-dim case
	if(dim == 2)
	{
		t_real angle = std::atan2(vec2[1], vec2[0])
			     - std::atan2(vec1[1], vec1[0]);

		return rotation_2d<t_mat>(angle);
	}

	// 3-dim and 4-dim homogeneous case
	else if(dim == 3 || (dim == 4 && force_3dim))
	{
		// get rotation axis from cross product
		t_vec axis = cross<t_vec>({ vec1, vec2 });
		t_real lenaxis = norm<t_vec>(axis);

		// rotation angle
		t_real angle = std::atan2(lenaxis, inner<t_vec>(vec1, vec2));

		// collinear vectors?
		if(equals<t_real>(angle, 0, eps))
			return unit<t_mat>(vec1.size());

		// antiparallel vectors?
		if(equals<t_real>(std::abs(angle), pi<t_real>, eps))
		{
			if(perp_vec)
			{
				axis = *perp_vec;
				lenaxis = norm<t_vec>(axis);
			}
			else
			{
				axis = perp<t_vec>(vec1, eps);
				lenaxis = t_real(1);
			}

			angle = pi<t_real>;
		}

		axis /= lenaxis;
		return rotation<t_mat, t_vec>(axis, angle, true);
	}

	// general case, equation (8) from (Zhelezov 2017)
	else
	{
		// matrix to rotate vec1 to [1, 0, 0, ...]
		t_mat mat1_x0 = rotation_x0<t_mat, t_vec, t_real>(vec1, eps);

		// matrix to rotate vec2 to [1, 0, 0, ...]
		t_mat mat2_x0 = rotation_x0<t_mat, t_vec, t_real>(vec2, eps);

		// matrix to rotate [1, 0, 0, ...] to vec2
		auto [mat2_x0_inv, mat2_ok] = inv<t_mat>(mat2_x0);

		// rotate vec1 to [1, 0, 0, ...], then rotate [1, 0, 0, ...] to vec2
		return mat2_x0_inv * mat1_x0;
	}
}


/**
 * extracts lines from polygon object, takes input from e.g. create_cube()
 * @returns [point pairs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
t_cont<t_vec> create_lines(const t_cont<t_vec>& vertices, const t_cont<t_cont<std::size_t>>& faces)
requires is_vec<t_vec>
{
	t_cont<t_vec> lineverts;

	auto line_already_seen = [&lineverts](const t_vec& vec1, const t_vec& vec2) -> bool
	{
		auto iter = lineverts.begin();

		while(1)
		{
			const t_vec& linevec1 = *iter;
			std::advance(iter, 1); if(iter == lineverts.end()) break;
			const t_vec& linevec2 = *iter;

			if(equals<t_vec>(vec1, linevec1) && equals<t_vec>(vec2, linevec2))
				return true;
			if(equals<t_vec>(vec1, linevec2) && equals<t_vec>(vec2, linevec1))
				return true;

			std::advance(iter, 1); if(iter == lineverts.end()) break;
		}

		return false;
	};

	for(const auto& face : faces)
	{
		// iterator to last point
		auto iter1 = face.begin();
		std::advance(iter1, face.size()-1);

		for(auto iter2 = face.begin(); iter2 != face.end(); std::advance(iter2, 1))
		{
			const t_vec& vec1 = vertices[*iter1];
			const t_vec& vec2 = vertices[*iter2];

			//if(!line_already_seen(vec1, vec2))
			{
				lineverts.push_back(vec1);
				lineverts.push_back(vec2);
			}

			iter1 = iter2;
		}
	}

	return lineverts;
}


/**
 * triangulates polygon object, takes input from e.g. create_cube()
 * @returns [triangles, face normals, vertex uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_vec>, t_cont<t_vec>>
create_triangles(const std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>,
	t_cont<t_vec>, t_cont<t_cont<t_vec>>>& tup)
requires is_vec<t_vec>
{
	const t_cont<t_vec>& vertices = std::get<0>(tup);
	const t_cont<t_cont<std::size_t>>& faces = std::get<1>(tup);
	const t_cont<t_vec>& normals = std::get<2>(tup);
	const t_cont<t_cont<t_vec>>& uvs = std::get<3>(tup);

	t_cont<t_vec> triangles;
	t_cont<t_vec> triag_normals;
	t_cont<t_vec> vert_uvs;

	auto iterFaces = faces.begin();
	auto iterNorms = normals.begin();
	auto iterUVs = uvs.begin();

	// iterate over faces
	while(iterFaces != faces.end())
	{
		// triangulate faces
		auto iterFaceVertIdx = iterFaces->begin();
		std::size_t vert1Idx = *iterFaceVertIdx;
		std::advance(iterFaceVertIdx, 1);
		std::size_t vert2Idx = *iterFaceVertIdx;

		const t_vec *puv1 = nullptr;
		const t_vec *puv2 = nullptr;
		const t_vec *puv3 = nullptr;

		typename t_cont<t_vec>::const_iterator iterFaceUVIdx;
		if(iterUVs != uvs.end() && iterFaceUVIdx != iterUVs->end())
		{
			iterFaceUVIdx = iterUVs->begin();

			puv1 = &(*iterFaceUVIdx);
			std::advance(iterFaceUVIdx, 1);
			puv2 = &(*iterFaceUVIdx);
		}

		// iterate over face vertices
		while(1)
		{
			std::advance(iterFaceVertIdx, 1);
			if(iterFaceVertIdx == iterFaces->end())
				break;
			std::size_t vert3Idx = *iterFaceVertIdx;

			if(iterUVs != uvs.end() && iterFaceUVIdx != iterUVs->end())
			{
				std::advance(iterFaceUVIdx, 1);
				puv3 = &(*iterFaceUVIdx);
			}

			// create triangle
			triangles.push_back(vertices[vert1Idx]);
			triangles.push_back(vertices[vert2Idx]);
			triangles.push_back(vertices[vert3Idx]);

			// triangle normal
			triag_normals.push_back(*iterNorms);
			//triag_normals.push_back(*iterNorms);
			//triag_normals.push_back(*iterNorms);

			// triangle vertex uvs
			if(puv1 && puv2 && puv3)
			{
				vert_uvs.push_back(*puv1);
				vert_uvs.push_back(*puv2);
				vert_uvs.push_back(*puv3);
			}


			// next vertex
			vert2Idx = vert3Idx;
			puv2 = puv3;
		}


		std::advance(iterFaces, 1);
		if(iterNorms != normals.end()) std::advance(iterNorms, 1);
		if(iterUVs != uvs.end()) std::advance(iterUVs, 1);
	}

	return std::make_tuple(triangles, triag_normals, vert_uvs);
}


/**
 * subdivides triangles
 * input: [triangle vertices, normals, uvs]
 * @returns [triangles, face normals, vertex uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_vec>, t_cont<t_vec>>
subdivide_triangles(const std::tuple<t_cont<t_vec>, t_cont<t_vec>, t_cont<t_vec>>& tup)
requires is_vec<t_vec>
{
	const t_cont<t_vec>& vertices = std::get<0>(tup);
	const t_cont<t_vec>& normals = std::get<1>(tup);
	const t_cont<t_vec>& uvs = std::get<2>(tup);

	t_cont<t_vec> vertices_new;
	t_cont<t_vec> normals_new;
	t_cont<t_vec> uvs_new;


	// iterate over triplets forming triangles
	auto itervert = vertices.begin();
	auto iternorm = normals.begin();
	auto iteruv = uvs.begin();

	while(itervert != vertices.end())
	{
		const t_vec& vec1 = *itervert;
		std::advance(itervert, 1); if(itervert == vertices.end()) break;
		const t_vec& vec2 = *itervert;
		std::advance(itervert, 1); if(itervert == vertices.end()) break;
		const t_vec& vec3 = *itervert;
		std::advance(itervert, 1);

		const t_vec vec12mid = avg<t_vec>({ vec1, vec2 });
		const t_vec vec23mid = avg<t_vec>({ vec2, vec3 });
		const t_vec vec31mid = avg<t_vec>({ vec3, vec1 });

		// triangle 1
		vertices_new.push_back(vec1);
		vertices_new.push_back(vec12mid);
		vertices_new.push_back(vec31mid);

		// triangle 2
		vertices_new.push_back(vec12mid);
		vertices_new.push_back(vec2);
		vertices_new.push_back(vec23mid);

		// triangle 3
		vertices_new.push_back(vec31mid);
		vertices_new.push_back(vec23mid);
		vertices_new.push_back(vec3);

		// triangle 4
		vertices_new.push_back(vec12mid);
		vertices_new.push_back(vec23mid);
		vertices_new.push_back(vec31mid);


		// duplicate normals for the four sub-triangles
		if(iternorm != normals.end())
		{
			normals_new.push_back(*iternorm);
			normals_new.push_back(*iternorm);
			normals_new.push_back(*iternorm);
			normals_new.push_back(*iternorm);

			std::advance(iternorm, 1);
		}


		// uv coords
		if(iteruv != uvs.end())
		{
			// uv coords at vertices
			const t_vec& uv1 = *iteruv;
			std::advance(iteruv, 1); if(iteruv == uvs.end()) break;
			const t_vec& uv2 = *iteruv;
			std::advance(iteruv, 1); if(iteruv == uvs.end()) break;
			const t_vec& uv3 = *iteruv;
			std::advance(iteruv, 1);

			const t_vec uv12mid = avg<t_vec>({ uv1, uv2 });
			const t_vec uv23mid = avg<t_vec>({ uv2, uv3 });
			const t_vec uv31mid = avg<t_vec>({ uv3, uv1 });

			// uvs of triangle 1
			uvs_new.push_back(uv1);
			uvs_new.push_back(uv12mid);
			uvs_new.push_back(uv31mid);

			// uvs of triangle 2
			uvs_new.push_back(uv12mid);
			uvs_new.push_back(uv2);
			uvs_new.push_back(uv23mid);

			// uvs of triangle 3
			uvs_new.push_back(uv31mid);
			uvs_new.push_back(uv23mid);
			uvs_new.push_back(uv3);

			// uvs of triangle 4
			uvs_new.push_back(uv12mid);
			uvs_new.push_back(uv23mid);
			uvs_new.push_back(uv31mid);
		}
	}

	return std::make_tuple(vertices_new, normals_new, uvs_new);
}


/**
 * subdivides triangles (with specified number of iterations)
 * input: [triangle vertices, normals, uvs]
 * @returns [triangles, face normals, vertex uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_vec>, t_cont<t_vec>>
subdivide_triangles(const std::tuple<t_cont<t_vec>, t_cont<t_vec>, t_cont<t_vec>>& tup, std::size_t iters)
requires is_vec<t_vec>
{
	auto tupDiv = tup;
	for(std::size_t i=0; i<iters; ++i)
		tupDiv = subdivide_triangles<t_vec, t_cont>(tupDiv);
	return tupDiv;
}


/**
 * create the faces of a sphere
 * input: [triangle vertices, normals, uvs] (like subdivide_triangles)
 * @returns [triangles, face normals, vertex uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_vec>, t_cont<t_vec>>
spherify(const std::tuple<t_cont<t_vec>, t_cont<t_vec>, t_cont<t_vec>>& tup,
	typename t_vec::value_type rad = 1)
requires is_vec<t_vec>
{
	const t_cont<t_vec>& vertices = std::get<0>(tup);
	//const t_cont<t_vec>& normals = std::get<1>(tup);
	const t_cont<t_vec>& uvs = std::get<2>(tup);

	t_cont<t_vec> vertices_new;
	t_cont<t_vec> normals_new;
	vertices_new.reserve(vertices.size());
	normals_new.reserve(vertices.size());


	// vertices
	for(t_vec vec : vertices)
	{
		vec /= norm<t_vec>(vec);
		vec *= rad;
		vertices_new.emplace_back(std::move(vec));
	}


	// face normals
	auto itervert = vertices.begin();
	// iterate over triplets forming triangles
	while(itervert != vertices.end())
	{
		const t_vec& vec1 = *itervert;
		std::advance(itervert, 1); if(itervert == vertices.end()) break;
		const t_vec& vec2 = *itervert;
		std::advance(itervert, 1); if(itervert == vertices.end()) break;
		const t_vec& vec3 = *itervert;
		std::advance(itervert, 1);

		t_vec vecmid = avg<t_vec>({ vec1, vec2, vec3 });
		vecmid /= norm<t_vec>(vecmid);
		normals_new.emplace_back(std::move(vecmid));
	}

	return std::make_tuple(vertices_new, normals_new, uvs);
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// 3-dim solids
// ----------------------------------------------------------------------------

/**
 * transforms vertices and normals using a matrix
 */
template<class t_mat, class t_vec, template<class...> class t_cont = std::vector>
void transform_obj(t_cont<t_vec>& verts, t_cont<t_vec>& norms, const t_mat& mat, bool is_3dhom=false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using size_t = decltype(mat.size1());

	// make sure a 3-vector and a 4-matrix are handled correctly in homogeneous coordinates
	if(is_3dhom && mat.size1()==4)
	{
		t_mat mat3 = mat;
		for(size_t i=0; i<3; ++i)
		{
			mat3(3,i) = 0;
			mat3(i,3) = 0;
		}
		mat3(3,3) = 1;

		for(auto& vert : verts)
		{
			vert = mat3 * vert;

			// add translation and normalise
			for(size_t i=0; i<3; ++i)
			{
				vert[i] += mat(i,3);
				vert[i] /= mat(3,3);
			}
		}

		for(auto& norm : norms)
			norm = mat3 * norm;
	}

	// standard case: just multiply
	else
	{
		for(auto& vert : verts)
			vert = mat * vert;

		for(auto& norm : norms)
			norm = mat * norm;
	}
}


/**
 * create a plane
 * @returns [vertices, face vertex indices, face normals, face uvs]
 */
template<class t_mat, class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_plane(const t_vec& norm, typename t_vec::value_type lx=1, typename t_vec::value_type ly=1)
requires is_vec<t_vec>
{
	t_vec norm_old = create<t_vec>({ 0, 0, -1 });
	t_vec rot_vec = create<t_vec>({ 1, 0, 0 });
	t_mat rot = rotation<t_mat, t_vec>(norm_old, norm, &rot_vec);

	t_cont<t_vec> vertices =
	{
		create<t_vec>({ -lx, -ly, 0. }),	// vertex 0
		create<t_vec>({ -lx, +ly, 0. }),	// vertex 1
		create<t_vec>({ +lx, +ly, 0. }),	// vertex 2
		create<t_vec>({ +lx, -ly, 0. }),	// vertex 3
	};

	// rotate according to given normal
	for(t_vec& vec : vertices)
		vec = rot * vec;

	t_cont<t_cont<std::size_t>> faces = { { 0, 1, 2, 3 } };
	t_cont<t_vec> normals = { norm };

	t_cont<t_cont<t_vec>> uvs =
	{{
		create<t_vec>({0, 0}),
		create<t_vec>({0, 1}),
		create<t_vec>({1, 1}),
		create<t_vec>({1, 0}),
	}};

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create a disk
 * @returns [vertices, face vertex indices, face normals, face uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_disk(typename t_vec::value_type r=1, std::size_t num_points=32)
requires is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	// vertices and uvs
	t_cont<t_vec> vertices;
	t_cont<t_vec> all_uvs;
	vertices.reserve(num_points);
	all_uvs.reserve(num_points);

	// vertices
	for(std::size_t pt=0; pt<num_points; ++pt)
	{
		const t_real phi = t_real(pt)/t_real(num_points) * t_real(2)*pi<t_real>;
		const t_real c = std::cos(phi);
		const t_real s = std::sin(phi);

		// outer vertices
		t_vec vert = create<t_vec>({ r*c, r*s, 0 });
		vertices.emplace_back(std::move(vert));

		t_vec uv = create<t_vec>({
			(t_real(1)+c) / t_real(2),
			(t_real(1)+s) / t_real(2) });
		all_uvs.emplace_back(std::move(uv));
	}

	// faces, normals & uvs
	t_cont<t_cont<std::size_t>> faces;
	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;

	t_cont<std::size_t> face(num_points);
	std::iota(face.begin(), face.end(), 0);
	faces.emplace_back(std::move(face));

	normals.emplace_back(create<t_vec>({0,0,1}));

	t_cont<t_vec> uv;
	uv.reserve(num_points);
	for(std::size_t vert_idx : faces[0])
		uv.push_back(all_uvs[vert_idx]);
	uvs.emplace_back(std::move(uv));

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create a cone
 * @returns [vertices, face vertex indices, face normals, face uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_cone(typename t_vec::value_type r = 1, typename t_vec::value_type h = 1,
	bool bWithCap = true, std::size_t num_points = 32)
requires is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	// vertices
	t_cont<t_vec> vertices;
	vertices.reserve(num_points + 1);

	// inner vertex
	vertices.emplace_back(create<t_vec>({ 0, 0, h }));

	for(std::size_t pt=0; pt<num_points; ++pt)
	{
		const t_real phi = t_real(pt)/t_real(num_points) * t_real(2)*pi<t_real>;
		const t_real c = std::cos(phi);
		const t_real s = std::sin(phi);

		// outer vertices
		t_vec vert = create<t_vec>({ r*c, r*s, 0 });
		vertices.emplace_back(std::move(vert));
	}

	// faces, normals & uvs
	t_cont<t_cont<std::size_t>> faces;
	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;	// TODO

	faces.reserve(num_points);
	normals.reserve(num_points);
	uvs.reserve(num_points);

	for(std::size_t face=0; face<num_points; ++face)
	{
		std::size_t idx0 = face + 1;	// outer 1
		std::size_t idx1 = (face == num_points-1 ? 1 : face + 2);	// outer 2
		std::size_t idx2 = 0;	// inner

		faces.push_back({ idx0, idx1, idx2 });

		t_vec n = cross<t_vec>({
			vertices[idx2]-vertices[idx0],
			vertices[idx1]-vertices[idx0] });
		n /= norm<t_vec>(n);

		normals.emplace_back(std::move(n));
	}

	if(bWithCap)
	{
		const auto [disk_vertices, disk_faces, disk_normals, disk_uvs]
			= create_disk<t_vec, t_cont>(r, num_points);

		// vertex indices have to be adapted for merging
		const std::size_t vert_start_idx = vertices.size();
		vertices.insert(vertices.end(), disk_vertices.begin(), disk_vertices.end());

		auto disk_faces_bottom = disk_faces;
		for(auto& disk_face : disk_faces_bottom)
		{
			for(auto& disk_face_idx : disk_face)
				disk_face_idx += vert_start_idx;
			std::reverse(disk_face.begin(), disk_face.end());
		}
		faces.insert(faces.end(), disk_faces_bottom.begin(), disk_faces_bottom.end());

		for(const auto& normal : disk_normals)
			normals.push_back(-normal);

		uvs.insert(uvs.end(), disk_uvs.begin(), disk_uvs.end());
	}

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create a cylinder
 * cyltype: 0 (no caps), 1 (with caps), 2 (arrow)
 * @returns [vertices, face vertex indices, face normals, face uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_cylinder(typename t_vec::value_type r = 1, typename t_vec::value_type h = 1,
	int cyltype = 0, std::size_t num_points = 32,
	typename t_vec::value_type arrow_r = 1.5, typename t_vec::value_type arrow_h = 0.5)
requires is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	// vertices
	t_cont<t_vec> vertices;
	t_cont<t_real> vertices_u;

	vertices.reserve(num_points*2);
	vertices_u.reserve(num_points*2);

	for(std::size_t pt=0; pt<num_points; ++pt)
	{
		const t_real u = t_real(pt)/t_real(num_points);
		const t_real phi = u * t_real(2)*pi<t_real>;
		const t_real c = std::cos(phi);
		const t_real s = std::sin(phi);

		t_vec top = create<t_vec>({ r*c, r*s, h*t_real(0.5) });
		t_vec bottom = create<t_vec>({ r*c, r*s, -h*t_real(0.5) });

		vertices.emplace_back(std::move(top));
		vertices.emplace_back(std::move(bottom));

		vertices_u.push_back(u);
	}

	// faces, normals & uvs
	t_cont<t_cont<std::size_t>> faces;
	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;

	faces.reserve(num_points);
	normals.reserve(num_points);
	uvs.reserve(num_points);

	for(std::size_t face=0; face<num_points; ++face)
	{
		std::size_t idx0 = face*2 + 0;	// top 1
		std::size_t idx1 = face*2 + 1;	// bottom 1
		std::size_t idx2 = (face >= num_points-1 ? 1 : face*2 + 3);	// bottom 2
		std::size_t idx3 = (face >= num_points-1 ? 0 : face*2 + 2);	// top 2

		t_vec n = cross<t_vec>({
			vertices[idx1]-vertices[idx0],
			vertices[idx3]-vertices[idx0] });
		n /= norm<t_vec>(n);

		faces.push_back({ idx0, idx1, idx2, idx3 });
		normals.emplace_back(std::move(n));

		t_real u1 = vertices_u[idx0];
		t_real u2 = (face >= num_points-1 ? 1 : vertices_u[idx3]);
		uvs.push_back({
			create<t_vec>({u1, 1}), create<t_vec>({u1, 0}),
			create<t_vec>({u2, 0}), create<t_vec>({u2, 1})
		});
	}

	if(cyltype > 0)
	{
		const auto [disk_vertices, disk_faces, disk_normals, disk_uvs]
			= create_disk<t_vec, t_cont>(r, num_points);

		// bottom lid
		// vertex indices have to be adapted for merging
		std::size_t vert_start_idx = vertices.size();
		const t_vec top = create<t_vec>({ 0, 0, h*t_real(0.5) });

		for(const auto& disk_vert : disk_vertices)
			vertices.push_back(disk_vert - top);

		auto disk_faces_bottom = disk_faces;
		for(auto& disk_face : disk_faces_bottom)
		{
			for(auto& disk_face_idx : disk_face)
				disk_face_idx += vert_start_idx;
			std::reverse(disk_face.begin(), disk_face.end());
		}
		faces.insert(faces.end(), disk_faces_bottom.begin(), disk_faces_bottom.end());

		for(const auto& normal : disk_normals)
			normals.push_back(-normal);

		uvs.insert(uvs.end(), disk_uvs.begin(), disk_uvs.end());


		vert_start_idx = vertices.size();

		if(cyltype == 1)	// top lid
		{
			for(const auto& disk_vert : disk_vertices)
				vertices.push_back(disk_vert + top);

			auto disk_faces_top = disk_faces;
			for(auto& disk_face : disk_faces_top)
				for(auto& disk_face_idx : disk_face)
					disk_face_idx += vert_start_idx;
			faces.insert(faces.end(), disk_faces_top.begin(), disk_faces_top.end());

			for(const auto& normal : disk_normals)
				normals.push_back(normal);

			uvs.insert(uvs.end(), disk_uvs.begin(), disk_uvs.end());
		}
		else if(cyltype == 2)	// arrow top
		{
			// no need to cap the arrow if the radii are equal
			bool bConeCap = !equals<t_real>(r, arrow_r);

			const auto [cone_vertices, cone_faces, cone_normals, cone_uvs] =
				create_cone<t_vec, t_cont>(arrow_r, arrow_h, bConeCap, num_points);

			for(const auto& cone_vert : cone_vertices)
				vertices.push_back(cone_vert + top);

			auto cone_faces_top = cone_faces;
			for(auto& cone_face : cone_faces_top)
				for(auto& cone_face_idx : cone_face)
					cone_face_idx += vert_start_idx;
			faces.insert(faces.end(), cone_faces_top.begin(), cone_faces_top.end());

			for(const auto& normal : cone_normals)
				normals.push_back(normal);

			uvs.insert(uvs.end(), cone_uvs.begin(), cone_uvs.end());
		}
	}

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create the faces of a cuboid
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_cuboid(typename t_vec::value_type lx=1, typename t_vec::value_type ly=1, typename t_vec::value_type lz=1)
requires is_vec<t_vec>
{
	t_cont<t_vec> vertices =
	{
		create<t_vec>({ +lx, -ly, -lz }),	// vertex 0
		create<t_vec>({ -lx, -ly, -lz }),	// vertex 1
		create<t_vec>({ -lx, +ly, -lz }),	// vertex 2
		create<t_vec>({ +lx, +ly, -lz }),	// vertex 3

		create<t_vec>({ -lx, -ly, +lz }),	// vertex 4
		create<t_vec>({ +lx, -ly, +lz }),	// vertex 5
		create<t_vec>({ +lx, +ly, +lz }),	// vertex 6
		create<t_vec>({ -lx, +ly, +lz }),	// vertex 7
	};

	t_cont<t_cont<std::size_t>> faces =
	{
		{ 0, 1, 2, 3 },	// -z face
		{ 4, 5, 6, 7 },	// +z face
		{ 1, 0, 5, 4 }, // -y face
		{ 7, 6, 3, 2 },	// +y face
		{ 1, 4, 7, 2 },	// -x face
		{ 5, 0, 3, 6 },	// +x face
	};

	t_cont<t_vec> normals =
	{
		create<t_vec>({ 0, 0, -1 }),	// -z face
		create<t_vec>({ 0, 0, +1 }),	// +z face
		create<t_vec>({ 0, -1, 0 }),	// -y face
		create<t_vec>({ 0, +1, 0 }),	// +y face
		create<t_vec>({ -1, 0, 0 }),	// -x face
		create<t_vec>({ +1, 0, 0 }),	// +x face
	};

	t_cont<t_cont<t_vec>> uvs =
	{
		// -z face
		{
			create<t_vec>({0, 0}),
			create<t_vec>({1, 0}),
			create<t_vec>({1, 1}),
			create<t_vec>({0, 1})
		},
		// +z face
		{
			create<t_vec>({0, 1}),
			create<t_vec>({1, 1}),
			create<t_vec>({1, 0}),
			create<t_vec>({0, 0})
		},
		// -y face
		{
			create<t_vec>({0, 1}),
			create<t_vec>({1, 1}),
			create<t_vec>({1, 0}),
			create<t_vec>({0, 0}),
		},
		// +y face
		{
			create<t_vec>({1, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 1}),
			create<t_vec>({1, 1})
		},
		// -x face
		{
			create<t_vec>({1, 1}),
			create<t_vec>({1, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 1})
		},
		// +x face
		{
			create<t_vec>({0, 0}),
			create<t_vec>({0, 1}),
			create<t_vec>({1, 1}),
			create<t_vec>({1, 0})
		},
	};

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create the faces of a cube
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_cube(typename t_vec::value_type l = 1)
requires is_vec<t_vec>
{
	return create_cuboid<t_vec, t_cont>(l, l, l);
}


/**
 * create the faces of a icosahedron
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_icosahedron(typename t_vec::value_type l = 1)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;
	const T g = golden<T>;

	t_cont<t_vec> vertices =
	{
		create<t_vec>({ 0, -l, -g*l }), create<t_vec>({ 0, -l, +g*l }),
		create<t_vec>({ 0, +l, -g*l }), create<t_vec>({ 0, +l, +g*l }),

		create<t_vec>({ -g*l, 0, -l }), create<t_vec>({ -g*l, 0, +l }),
		create<t_vec>({ +g*l, 0, -l }), create<t_vec>({ +g*l, 0, +l }),

		create<t_vec>({ -l, -g*l, 0 }), create<t_vec>({ -l, +g*l, 0 }),
		create<t_vec>({ +l, -g*l, 0 }), create<t_vec>({ +l, +g*l, 0 }),
	};

	t_cont<t_cont<std::size_t>> faces =
	{
		{ 4, 2, 0 }, { 0, 6, 10 }, { 10, 7, 1 }, { 1, 3, 5 }, { 5, 9, 4 },
		{ 7, 10, 6 }, { 6, 0, 2 }, { 2, 4, 9 }, { 9, 5, 3 }, { 3, 1, 7 },
		{ 0, 10, 8 }, { 10, 1, 8 }, { 1, 5, 8 }, { 5, 4, 8 }, { 4, 0, 8 },
		{ 3, 7, 11 }, { 7, 6, 11 }, { 6, 2, 11 }, { 2, 9, 11 }, { 9, 3, 11 },
	};


	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;

	normals.reserve(faces.size());
	uvs.reserve(faces.size());

	for(const auto& face : faces)
	{
		auto iter = face.begin();
		const t_vec& vec1 = *(vertices.begin() + *iter); std::advance(iter,1);
		const t_vec& vec2 = *(vertices.begin() + *iter); std::advance(iter,1);
		const t_vec& vec3 = *(vertices.begin() + *iter);

		const t_vec vec12 = vec2 - vec1;
		const t_vec vec13 = vec3 - vec1;

		// TODO
		uvs.emplace_back(t_cont<t_vec>{{
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}) }});

		t_vec n = cross<t_vec>({vec12, vec13});
		n /= norm<t_vec>(n);
		normals.emplace_back(std::move(n));
	}

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create the faces of a dodecahedron
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_dodecahedron(typename t_vec::value_type l = 1)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;
	const T g = golden<T>;

	t_cont<t_vec> vertices =
	{
		create<t_vec>({ l, l, l }), create<t_vec>({ l, l, -l }),
		create<t_vec>({ l, -l, l }), create<t_vec>({ l, -l, -l }),

		create<t_vec>({ -l, l, l }), create<t_vec>({ -l, l, -l }),
		create<t_vec>({ -l, -l, l }), create<t_vec>({ -l, -l, -l }),

		create<t_vec>({ 0, T{l}/g, g*l }), create<t_vec>({ 0, T{l}/g, -g*l }),
		create<t_vec>({ 0, -T{l}/g, g*l }), create<t_vec>({ 0, -T{l}/g, -g*l }),

		create<t_vec>({ g*l, 0, T{l}/g }), create<t_vec>({ g*l, 0, -T{l}/g }),
		create<t_vec>({ -g*l, 0, T{l}/g }), create<t_vec>({ -g*l, 0, -T{l}/g }),

		create<t_vec>({ T{l}/g, g*l, 0 }), create<t_vec>({ T{l}/g, -g*l, 0 }),
		create<t_vec>({ -T{l}/g, g*l, 0 }), create<t_vec>({ -T{l}/g, -g*l, 0 }),
	};

	t_cont<t_cont<std::size_t>> faces =
	{
		{ 0, 16, 18, 4, 8 }, { 0, 8, 10, 2, 12 }, { 0, 12, 13, 1, 16 },
		{ 1, 9, 5, 18, 16 }, { 1, 13, 3, 11, 9 }, { 2, 17, 3, 13, 12 },
		{ 3, 17, 19, 7, 11 }, { 2, 10, 6, 19, 17 }, { 4, 14, 6, 10, 8 },
		{ 4, 18, 5, 15, 14 }, { 5, 9, 11, 7, 15 }, { 6, 14, 15, 7, 19 },
	};


	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;

	normals.reserve(faces.size());
	uvs.reserve(faces.size());

	for(const auto& face : faces)
	{
		auto iter = face.begin();
		const t_vec& vec1 = *(vertices.begin() + *iter); std::advance(iter,1);
		const t_vec& vec2 = *(vertices.begin() + *iter); std::advance(iter,1);
		const t_vec& vec3 = *(vertices.begin() + *iter);

		const t_vec vec12 = vec2 - vec1;
		const t_vec vec13 = vec3 - vec1;

		// TODO
		uvs.emplace_back(t_cont<t_vec>{{
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}) }});

		t_vec n = cross<t_vec>({vec12, vec13});
		n /= norm<t_vec>(n);
		normals.emplace_back(std::move(n));
	}

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create the faces of a octahedron
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_octahedron(typename t_vec::value_type l = 1)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_cont<t_vec> vertices =
	{
		create<t_vec>({ +l, 0, 0 }),	// vertex 0
		create<t_vec>({ 0, +l, 0 }),	// vertex 1
		create<t_vec>({ 0, 0, +l }),	// vertex 2

		create<t_vec>({ -l, 0, 0 }),	// vertex 3
		create<t_vec>({ 0, -l, 0 }),	// vertex 4
		create<t_vec>({ 0, 0, -l }),	// vertex 5
	};

	t_cont<t_cont<std::size_t>> faces =
	{
		{ 2, 0, 1 }, { 0, 5, 1 }, { 5, 3, 1 }, { 3, 2, 1 },	// upper half
		{ 0, 2, 4 }, { 5, 0, 4 }, { 3, 5, 4 }, { 2, 3, 4 },	// lower half
	};


	const T len = std::sqrt(3);

	t_cont<t_vec> normals =
	{
		create<t_vec>({ +1/len, +1/len, +1/len }),
		create<t_vec>({ +1/len, +1/len, -1/len }),
		create<t_vec>({ -1/len, +1/len, -1/len }),
		create<t_vec>({ -1/len, +1/len, +1/len }),

		create<t_vec>({ +1/len, -1/len, +1/len }),
		create<t_vec>({ +1/len, -1/len, -1/len }),
		create<t_vec>({ -1/len, -1/len, -1/len }),
		create<t_vec>({ -1/len, -1/len, +1/len }),
	};

	t_cont<t_cont<t_vec>> uvs =
	{
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },

		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
	};

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create the faces of a tetrahedron
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_tetrahedron(typename t_vec::value_type l = 1)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_cont<t_vec> vertices =
	{
		create<t_vec>({ -l, -l, +l }),	// vertex 0
		create<t_vec>({ +l, +l, +l }),	// vertex 1
		create<t_vec>({ -l, +l, -l }),	// vertex 2
		create<t_vec>({ +l, -l, -l }),	// vertex 3
	};

	t_cont<t_cont<std::size_t>> faces =
	{
		{ 1, 2, 0 }, { 2, 1, 3 },	// connected to upper edge
		{ 0, 3, 1 }, { 3, 0, 2 },	// connected to lower edge
	};


	const T len = std::sqrt(3);

	t_cont<t_vec> normals =
	{
		create<t_vec>({ -1/len, +1/len, +1/len }),
		create<t_vec>({ +1/len, +1/len, -1/len }),
		create<t_vec>({ +1/len, -1/len, +1/len }),
		create<t_vec>({ -1/len, -1/len, -1/len }),
	};

	t_cont<t_cont<t_vec>> uvs =
	{
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
	};

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * calculates the bounding sphere of a collection of vertices
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_vec, typename t_vec::value_type> bounding_sphere(
	const t_cont<t_vec>& verts)
requires is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	t_real rad{};
	t_vec center = mean<t_vec, t_cont>(verts);

	for(const t_vec& vec : verts)
	{
		t_vec vecCur = vec-center;

		t_real dot = inner<t_vec>(vecCur, vecCur);
		rad = std::max(rad, dot);
	}

	rad = std::sqrt(rad);
	return std::make_tuple(center, rad);
}


/**
 * calculates the bounding box of a collection of vertices
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_vec, t_vec> bounding_box(
	const t_cont<t_vec>& verts,
	const t_vec* oldmin = nullptr, const t_vec* oldmax = nullptr)
requires is_vec<t_vec>
{
	if(verts.size() == 0)
		return std::make_tuple(t_vec{}, t_vec{});

	using t_real = typename t_vec::value_type;
	const std::size_t dim = verts[0].size();
	if(dim == 0)
		return std::make_tuple(t_vec{}, t_vec{});

	t_vec min, max;
	if(oldmin && oldmax)
	{
		min = *oldmin;
		max = *oldmax;
	}
	else
	{
		min = create<t_vec>(dim);
		max = create<t_vec>(dim);

		for(std::size_t i=0; i<dim; ++i)
		{
			min[i] = std::numeric_limits<t_real>::max();
			max[i] = std::numeric_limits<t_real>::lowest();
		}
	}

	for(const t_vec& vert : verts)
	{
		for(std::size_t i=0; i<dim; ++i)
		{
			min[i] = std::min(min[i], vert[i]);
			max[i] = std::max(max[i], vert[i]);
		}
	}

	return std::make_tuple(min, max);
}


/**
 * calculates the bounding box of a collection of vertices
 * @see https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_vec, t_vec> bounding_box(
	const t_cont<t_cont<t_vec>>& allverts, std::size_t dim = 3)
requires is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	t_vec min = create<t_vec>(dim);
	t_vec max = create<t_vec>(dim);

	for(std::size_t i=0; i<dim; ++i)
	{
		min[i] = std::numeric_limits<t_real>::max();
		max[i] = std::numeric_limits<t_real>::lowest();
	}

	for(const t_cont<t_vec>& verts : allverts)
		std::tie(min, max) = bounding_box<t_vec, t_cont>(verts, &min, &max);

	return std::make_tuple(min, max);
}


/**
 * calculates the bounding box of a sphere
 */
template<class t_vec, class t_real=typename t_vec::value_type> requires is_vec<t_vec>
std::tuple<t_vec, t_vec> sphere_bounding_box(const t_vec& pos, t_real rad)
{
	const std::size_t dim = pos.size();
	if(dim == 0)
		return std::make_tuple(t_vec{}, t_vec{});

	t_vec min = create<t_vec>(dim);
	t_vec max = create<t_vec>(dim);

	for(std::size_t i=0; i<dim; ++i)
	{
		min[i] = pos[i] - rad;
		max[i] = pos[i] + rad;
	}

	return std::make_tuple(min, max);
}


/**
 * calculates the bounding box of a collection of spheres
 */
template<class t_vec, template<class...> class t_cont = std::vector,
	class t_real = typename t_vec::value_type>
std::tuple<t_vec, t_vec> sphere_bounding_box(
	const t_cont<std::tuple<t_vec, t_real>>& spheres, std::size_t dim = 3)
requires is_vec<t_vec>
{
	t_vec min = create<t_vec>(dim);
	for(std::size_t i=0; i<dim; ++i)
		min[i] = std::numeric_limits<t_real>::max();
	t_vec max = -min;

	for(const auto& sphere : spheres)
	{
		auto [spheremin, spheremax] = sphere_bounding_box<t_vec, t_real>(
			std::get<0>(sphere), std::get<1>(sphere));

		for(std::size_t i=0; i<dim; ++i)
		{
			min[i] = std::min(min[i], spheremin[i]);
			max[i] = std::max(max[i], spheremax[i]);
		}
	}

	return std::make_tuple(min, max);
}


/**
 * checks if a point is inside a bounding box
 * @see https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection
 */
template<class t_vec> requires is_vec<t_vec>
bool in_bounding_box(
	const t_vec& pt, const std::tuple<t_vec, t_vec>& bb)
{
	const std::size_t dim = pt.size();

	const t_vec& min = std::get<0>(bb);
	const t_vec& max = std::get<1>(bb);

	for(std::size_t i=0; i<dim; ++i)
	{
		if(pt[i] < min[i] || pt[i] > max[i])
			return false;
	}

	return true;
}


/**
 * checks for bounding box intersection
 * @see https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection
 */
template<class t_vec> requires is_vec<t_vec>
bool collide_bounding_boxes(
	const std::tuple<t_vec, t_vec>& bb1,
	const std::tuple<t_vec, t_vec>& bb2)
{
	// invalid bounding boxes?
	if(std::get<0>(bb1).size()==0 || std::get<1>(bb1).size()==0 ||
		std::get<0>(bb2).size()==0 || std::get<1>(bb2).size()==0)
		return false;

	const std::size_t dim = std::min(std::get<0>(bb1).size(), std::get<0>(bb2).size());
	for(std::size_t i=0; i<dim; ++i)
	{
		if(std::get<0>(bb1)[i] > std::get<1>(bb2)[i])
			return false;
		if(std::get<0>(bb2)[i] > std::get<1>(bb1)[i])
			return false;
	}
	return true;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// 3-dim algos in homogeneous coordinates
// ----------------------------------------------------------------------------

/**
 * does the homogeneous matrix mat have a translation component?
 */
template<class t_mat, class t_real = typename t_mat::value_type>
bool hom_has_translation_components(const t_mat& mat,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat>
{
	const std::size_t N = mat.size1();
	if(N != mat.size2())
		return false;

	// translation?
	for(std::size_t i=0; i<N-1; ++i)
	{
		if(!equals<t_real>(mat(i, N-1), t_real(0), eps))
			return true;
	}

	return false;
}


/**
 * is mat a centering matrix in homogeneous coords?
 */
template<class t_mat, class t_real = typename t_mat::value_type>
bool hom_is_centring(const t_mat& mat,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat>
{
	const std::size_t N = mat.size1();
	if(N != mat.size2())
		return false;

	// is the left-upper 3x3 a rotation matrix (and no unit matrix)?
	if(!is_unit(submat<t_mat>(mat, 0, 0, 3, 3), eps))
		return false;

	// translation?
	if(hom_has_translation_components<t_mat, t_real>(mat, eps))
		return true;

	return false;
}


/**
 * project a homogeneous vector to screen coordinates
 * @returns [vecPersp, vecScreen]
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluProject.xml
 */
template<class t_mat, class t_vec>
std::tuple<t_vec, t_vec> hom_to_screen_coords(const t_vec& vec4,
	const t_mat& matModelView, const t_mat& matProj, const t_mat& matViewport,
	bool bFlipY = false, bool bFlipX = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	// perspective trafo and divide
	t_vec vecPersp = matProj * matModelView * vec4;
	vecPersp /= vecPersp[3];

	// viewport trafo
	t_vec vec = matViewport * vecPersp;

	// flip y coordinate
	if(bFlipY) vec[1] = matViewport(1,1)*2 - vec[1];
	// flip x coordinate
	if(bFlipX) vec[0] = matViewport(0,0)*2 - vec[0];

	return std::make_tuple(vecPersp, vec);
}


/**
 * calculate world coordinates from screen coordinates
 * (vary zPlane to get the points of the z-line at constant (x,y))
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluUnProject.xml
 */
template<class t_mat, class t_vec>
t_vec hom_from_screen_coords(
	typename t_vec::value_type xScreen, typename t_vec::value_type yScreen, typename t_vec::value_type zPlane,
	const t_mat& matModelView_inv, const t_mat& matProj_inv, const t_mat& matViewport_inv,
	const t_mat* pmatViewport = nullptr, bool bFlipY = false, bool bFlipX = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_vec vecScreen = create<t_vec>({xScreen, yScreen, zPlane, 1.});

	// flip y coordinate
	if(pmatViewport && bFlipY) vecScreen[1] = (*pmatViewport)(1,1)*2 - vecScreen[1];
	// flip x coordinate
	if(pmatViewport && bFlipX) vecScreen[0] = (*pmatViewport)(0,0)*2 - vecScreen[0];

	t_vec vecWorld = matModelView_inv * matProj_inv * matViewport_inv * vecScreen;

	vecWorld /= vecWorld[3];
	return vecWorld;
}


/**
 * calculate line from screen coordinates
 * @returns [pos, dir]
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluUnProject.xml
 */
template<class t_mat, class t_vec>
std::tuple<t_vec, t_vec> hom_line_from_screen_coords(
	typename t_vec::value_type xScreen, typename t_vec::value_type yScreen,
	typename t_vec::value_type z1, typename t_vec::value_type z2,
	const t_mat& matModelView_inv, const t_mat& matProj_inv, const t_mat& matViewport_inv,
	const t_mat* pmatViewport = nullptr, bool bFlipY = false, bool bFlipX = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	const t_vec lineOrg = hom_from_screen_coords<t_mat, t_vec>(
		xScreen, yScreen, z1, matModelView_inv, matProj_inv,
		matViewport_inv, pmatViewport, bFlipY, bFlipX);
	const t_vec linePos2 = hom_from_screen_coords<t_mat, t_vec>(
		xScreen, yScreen, z2, matModelView_inv, matProj_inv,
		matViewport_inv, pmatViewport, bFlipY, bFlipX);

	t_vec lineDir = linePos2 - lineOrg;
	lineDir /= norm<t_vec>(lineDir);

	return std::make_tuple(lineOrg, lineDir);
}


/**
 * perspective matrix (homogeneous 4x4)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluPerspective.xml
 */
template<class t_mat, typename t_scalar=typename t_mat::value_type>
t_mat hom_perspective(
	t_scalar n = 0.01, t_scalar f = 100.,
	t_scalar fov = 0.5*pi<typename t_mat::value_type>, t_scalar ratio = 3./4.,
	bool bLHS = true, bool bZ01 = false)
requires is_mat<t_mat>
{
	const t_scalar c = 1./std::tan(0.5 * fov);
	const t_scalar n0 = bZ01 ? t_scalar{0} : n;
	const t_scalar sc = bZ01 ? t_scalar{1} : t_scalar{2};
	const t_scalar zs = bLHS ? t_scalar{-1} : t_scalar{1};
	const t_scalar range_nf = std::abs(f-n);

	//         ( x*c*r                           )      ( -x*c*r/z                         )
	//         ( y*c                             )      ( -y*c/z                           )
	// P * x = ( z*(n0+f)/(n-f) + w*sc*n*f/(n-f) )  =>  ( -(n0+f)/(n-f) - w/z*sc*n*f/(n-f) )
	//         ( -z                              )      ( 1                                )
	return create<t_mat>({
		c*ratio,  0.,  0.,                   0.,
		0,        c,   0.,                   0.,
		0.,       0.,  zs*(n0+f)/range_nf,   -sc*n*f/range_nf,
		0.,       0.,  zs,                   0.
	});
}


/**
 * orthographic projection matrix (homogeneous 4x4)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glOrtho.xml
 * @see https://en.wikipedia.org/wiki/Orthographic_projection
 */
template<class t_mat, typename t_scalar=typename t_mat::value_type>
t_mat hom_ortho_sym(
	t_scalar n = 0.01, t_scalar f = 100.,
	t_scalar range_lr = 2., t_scalar range_bt = 2.,
	bool bLHS = true, bool bMap05 = false)
requires is_mat<t_mat>
{
	// map ranges into [-0.5, 0.5] or [-1, 1] else
	const t_scalar sc = bMap05 ? t_scalar{1} : t_scalar{2};
	const t_scalar zs = bLHS ? t_scalar{-1} : t_scalar{1};

	const t_scalar range_nf = std::abs(f) + std::abs(n);

	// centring
	const t_scalar tr_z = sc*t_scalar{0.5} * (n+f) / range_nf;

	//         ( sc_x*x        )
	//         ( sc_y*y        )
	// P * x = ( sc_z*z - tr_z )
	//         ( 1             )
	return create<t_mat>({
		sc/range_lr, 0.,          0.,              0.,
		0.,          sc/range_bt, 0.,              0.,
		0.,          0.,          zs*sc/range_nf,  zs*tr_z,
		0.,          0.,          0.,              1.
	});
}


/**
 * orthographic projection matrix (homogeneous 4x4)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glOrtho.xml
 * @see https://en.wikipedia.org/wiki/Orthographic_projection
 */
template<class t_mat, typename t_scalar=typename t_mat::value_type>
t_mat hom_ortho(
	t_scalar n = 0.01, t_scalar f = 100.,
	t_scalar l = -1., t_scalar r = 1.,
	t_scalar b = -1., t_scalar t = 1.,
	bool bLHS = true, bool bMap05 = false)
requires is_mat<t_mat>
{
	// map ranges into [-0.5, 0.5] or [-1, 1] else
	const t_scalar sc = bMap05 ? t_scalar{1} : t_scalar{2};
	const t_scalar zs = bLHS ? t_scalar{-1} : t_scalar{1};

	const t_scalar range_lr = r - l;
	const t_scalar range_bt = t - b;
	const t_scalar range_nf = f - n;

	// centring
	const t_scalar tr_x = sc*t_scalar{0.5} * (l+r) / range_lr;
	const t_scalar tr_y = sc*t_scalar{0.5} * (b+t) / range_bt;
	const t_scalar tr_z = sc*t_scalar{0.5} * (n+f) / range_nf;

	//         ( sc_x*x - tr_x )
	//         ( sc_y*y - tr_y )
	// P * x = ( sc_z*z - tr_z )
	//         ( 1             )
	return create<t_mat>({
		sc/range_lr, 0.,          0.,             -tr_x,
		0.,          sc/range_bt, 0.,             -tr_y,
		0.,          0.,          zs*sc/range_nf, zs*tr_z,
		0.,          0.,          0.,             1.
	});
}


/**
 * viewport matrix (homogeneous 4x4)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glViewport.xml
 */
template<class t_mat, typename t_scalar=typename t_mat::value_type>
t_mat hom_viewport(t_scalar w, t_scalar h, t_scalar n = 0, t_scalar f = 1)
requires is_mat<t_mat>
{
	t_scalar sc{0.5};

	return create<t_mat>({
		sc*w,  0.,    0.,        sc*w,
		0,     sc*h,  0.,        sc*h,
		0.,    0.,    sc*(f-n),  sc*(f+n),
		0.,    0.,    0.,        1.
	});
}


/**
 * translation matrix in homogeneous coordinates
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat hom_translation(t_real x, t_real y, t_real z)
requires is_mat<t_mat>
{
	return create<t_mat>({
		1.,  0.,  0.,  x,
		0.,  1.,  0.,  y,
		0.,  0.,  1.,  z,
		0.,  0.,  0.,  1.
	});
}


/**
 * scaling matrix in homogeneous coordinates
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat hom_scaling(t_real x, t_real y, t_real z)
requires is_mat<t_mat>
{
	return create<t_mat>({
		x,   0.,  0.,  0.,
		0.,  y,   0.,  0.,
		0.,  0.,  z,   0.,
		0.,  0.,  0.,  1.
	});
}


/**
 * mirror matrix in homogeneous coordinates
 */
template<class t_mat, class t_vec>
t_mat hom_mirror(const t_vec& axis, bool is_normalised=1)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat mat = ortho_mirror_op<t_mat, t_vec>(axis, is_normalised);

	t_mat mat_hom = unit<t_mat>(4,4);
	tl2::convert<t_mat, t_mat>(mat_hom, mat);

	// in case the matrix is statically sized and was already larger than 4x4
	mat_hom(3, 3) = 1.;
	return mat_hom;
}


/**
 * mirror matrix in homogeneous coordinates with translation
 */
template<class t_mat, class t_vec>
t_mat hom_mirror(const t_vec& axis, const t_vec& pos, bool is_normalised=1)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat mirr = hom_mirror<t_mat, t_vec>(axis, is_normalised);

	t_mat offs = hom_translation<t_mat, t_vec>(pos);
	t_mat offs_inv = hom_translation<t_mat, t_vec>(-pos);

	return offs * mirr * offs_inv;
}


/**
 * "look at" matrix in homogeneous coordinates
 * @see (Sellers 2014), pp. 78-79
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluLookAt.xml
 */
template<class t_mat, class t_vec>
t_mat hom_lookat(const t_vec& pos, const t_vec& target, const t_vec& _up)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_real = typename t_mat::value_type;

	// create orthonormal system
	t_vec dir = -(target - pos);
	dir = dir / norm<t_vec>(dir);

	t_vec side = cross<t_vec>({_up, dir});
	side = side / norm<t_vec>(side);

	t_vec up = cross<t_vec>({dir, side});
	//up = up / norm<t_vec>(up);

	// inverted/transposed rotation matrix
	t_mat rot_inv = unit<t_mat>(4);
	set_row<t_mat, t_vec>(rot_inv, side, 0);
	set_row<t_mat, t_vec>(rot_inv, up, 1);
	set_row<t_mat, t_vec>(rot_inv, dir, 2);

	// inverted translation matrix
	t_mat trans_inv = hom_translation<t_mat, t_real>(
		-pos[0], -pos[1], -pos[2]);

	return rot_inv * trans_inv;
}


/**
 * shear matrix
 * @see https://en.wikipedia.org/wiki/Shear_matrix
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat shear(std::size_t ROWS, std::size_t COLS, std::size_t i, std::size_t j, t_real s)
requires is_mat<t_mat>
{
	t_mat mat = unit<t_mat>(ROWS, COLS);
	mat(i,j) = s;
	return mat;
}


/**
 * rotation matrix in homogeneous coordinates
 */
template<class t_mat, class t_vec>
t_mat hom_rotation(const t_vec& vec1, const t_vec& vec2,
	const t_vec *perp_vec = nullptr)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat rot = rotation<t_mat, t_vec>(vec1, vec2, perp_vec);

	return create<t_mat>({
		rot(0,0), rot(0,1), rot(0,2), 0.,
		rot(1,0), rot(1,1), rot(1,2), 0.,
		rot(2,0), rot(2,1), rot(2,2), 0.,
		0.,       0.,       0.,       1.
	});
}


/**
 * rotation matrix in homogeneous coordinates
 */
template<class t_mat, class t_vec>
t_mat hom_rotation(const t_vec& axis, typename t_vec::value_type angle,
	bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat rot = rotation<t_mat, t_vec>(axis, angle, is_normalised);

	return create<t_mat>({
		rot(0,0), rot(0,1), rot(0,2), 0.,
		rot(1,0), rot(1,1), rot(1,2), 0.,
		rot(2,0), rot(2,1), rot(2,2), 0.,
		0.,       0.,       0.,       1.
	});
}

// ----------------------------------------------------------------------------



/**
 * arrow matrix
 */
template<class t_vec, class t_mat, class t_real = typename t_vec::value_type>
t_mat get_arrow_matrix(
	const t_vec& vecTo,
	t_real postscale = 1, const t_vec& vecPostTrans = create<t_vec>({0,0,0.5}),
	const t_vec& vecFrom = create<t_vec>({0,0,1}),
	t_real prescale = 1, const t_vec& vecPreTrans = create<t_vec>({0,0,0}),
	const t_vec* perp_vec = nullptr)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat mat = unit<t_mat>(4);

	mat *= hom_translation<t_mat>(vecPreTrans[0], vecPreTrans[1], vecPreTrans[2]);
	mat *= hom_scaling<t_mat>(prescale, prescale, prescale);

	mat *= hom_rotation<t_mat, t_vec>(vecFrom, vecTo, perp_vec);

	mat *= hom_scaling<t_mat>(postscale, postscale, postscale);
	mat *= hom_translation<t_mat>(vecPostTrans[0], vecPostTrans[1], vecPostTrans[2]);

	return mat;
}



// ----------------------------------------------------------------------------
// complex algos
// ----------------------------------------------------------------------------

/**
 * split a complex vector into two vectors with the real and imag parts
 */
template<class t_vec_cplx, class t_vec_real>
std::tuple<t_vec_real, t_vec_real> split_cplx(const t_vec_cplx& vec)
requires is_complex<typename t_vec_cplx::value_type> && is_vec<t_vec_cplx> && is_vec<t_vec_real>
{
	using t_real = typename t_vec_real::value_type;

	t_vec_real vecRe = zero<t_vec_real>(vec.size());
	t_vec_real vecIm = zero<t_vec_real>(vec.size());

	auto iter = vec.begin();
	auto iterRe = vecRe.begin();
	auto iterIm = vecIm.begin();

	for(; iter!=vec.end(); )
	{
		*iterRe = t_real{iter->real()};
		*iterIm = t_real{iter->imag()};

		std::advance(iter, 1);
		std::advance(iterRe, 1);
		std::advance(iterIm, 1);
	}

	return std::make_tuple(vecRe, vecIm);
}


/**
 * split a complex matrix into two matrices with the real and imag parts
 */
template<class t_mat_cplx, class t_mat_real>
std::tuple<t_mat_real, t_mat_real> split_cplx(const t_mat_cplx& mat)
requires is_complex<typename t_mat_cplx::value_type> && is_mat<t_mat_cplx> && is_mat<t_mat_real>
{
	t_mat_real matRe = zero<t_mat_real>(mat.size1(), mat.size2());
	t_mat_real matIm = zero<t_mat_real>(mat.size1(), mat.size2());

	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			matRe(i,j) = mat(i,j).real();
			matIm(i,j) = mat(i,j).imag();
		}
	}

	return std::make_tuple(matRe, matIm);
}


/**
 * split a complex vector into a real vector with the real and imag components
 */
template<class t_vec_cplx, class t_vec_real>
t_vec_real split_cplx_real(const t_vec_cplx& vec)
requires is_complex<typename t_vec_cplx::value_type> && is_vec<t_vec_cplx> && is_vec<t_vec_real>
{
	using t_real = typename t_vec_real::value_type;
	t_vec_real vecReal = zero<t_vec_real>(vec.size() * 2);

	auto iter = vec.begin();
	auto iterRe = vecReal.begin();
	auto iterIm = std::next(vecReal.begin(), vec.size());

	for(; iter!=vec.end(); )
	{
		*iterRe = t_real{iter->real()};
		*iterIm = t_real{iter->imag()};

		std::advance(iter, 1);
		std::advance(iterRe, 1);
		std::advance(iterIm, 1);
	}

	return vecReal;
}


/**
 * split a complex matrix into a real matrix with the real and imag components
 */
template<class t_mat_cplx, class t_mat_real>
t_mat_real split_cplx_real(const t_mat_cplx& mat)
requires is_complex<typename t_mat_cplx::value_type> && is_mat<t_mat_cplx> && is_mat<t_mat_real>
{
	t_mat_real matRe = zero<t_mat_real>(mat.size1()*2, mat.size2()*2);

	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			matRe(i,j) = mat(i,j).real();
			matRe(i + mat.size1(), j + mat.size2()) = mat(i,j).imag();
		}
	}

	return matRe;
}


/**
 * create a complex vector from two real vectors
 */
template<class t_vec_cplx, class t_vec_real>
t_vec_cplx unite_cplx(const t_vec_real& vecRe, const t_vec_real& vecIm)
requires is_complex<typename t_vec_cplx::value_type> && is_vec<t_vec_cplx> && is_vec<t_vec_real>
{
	assert(vecRe.size() == vecIm.size());
	t_vec_cplx vec = zero<t_vec_cplx>(vecRe.size());

	auto iter = vec.begin();
	auto iterRe = vecRe.begin();
	auto iterIm = vecIm.begin();

	for(; iterRe!=vecRe.end(); )
	{
		(*iter).real(*iterRe);
		(*iter).imag(*iterIm);

		std::advance(iter, 1);
		std::advance(iterRe, 1);
		std::advance(iterIm, 1);
	}

	return vec;
}


/**
 * create a complex matrix from two real matrices
 */
template<class t_mat_cplx, class t_mat_real>
t_mat_cplx unite_cplx(const t_mat_real& matRe, const t_mat_real& matIm)
requires is_complex<typename t_mat_cplx::value_type> && is_mat<t_mat_cplx> && is_mat<t_mat_real>
{
	assert(matRe.size1() == matIm.size1());
	assert(matRe.size2() == matIm.size2());

	t_mat_cplx mat = zero<t_mat_cplx>(matRe.size1(), matRe.size2());

	for(std::size_t i=0; i<matRe.size1(); ++i)
	{
		for(std::size_t j=0; j<matRe.size2(); ++j)
		{
			mat(i, j).real(matRe(i, j));
			mat(i, j).imag(matIm(i, j));
		}
	}

	return mat;
}


/**
 * create a complex vector from a real vectors having real and imag components
 */
template<class t_vec_cplx, class t_vec_real>
t_vec_cplx unite_cplx_real(const t_vec_real& vecReIm)
requires is_complex<typename t_vec_cplx::value_type> && is_vec<t_vec_cplx> && is_vec<t_vec_real>
{
	t_vec_cplx vec = zero<t_vec_cplx>(vecReIm.size() / 2);

	auto iter = vec.begin();
	auto iterRe = vecReIm.begin();
	auto iterIm = std::next(vecReIm.begin(), vec.size());

	for(; iter!=vec.end(); )
	{
		(*iter).real(*iterRe);
		(*iter).imag(*iterIm);

		std::advance(iter, 1);
		std::advance(iterRe, 1);
		std::advance(iterIm, 1);
	}

	return vec;
}


/**
 * create a complex matrix from a real matrix having real and imag components
 */
template<class t_mat_cplx, class t_mat_real>
t_mat_cplx unite_cplx_real(const t_mat_real& matReIm)
requires is_complex<typename t_mat_cplx::value_type> && is_mat<t_mat_cplx> && is_mat<t_mat_real>
{
	t_mat_cplx mat = zero<t_mat_cplx>(matReIm.size1()/2, matReIm.size2()/2);

	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			mat(i, j).real(matReIm(i, j));
			mat(i, j).imag(matReIm(i+mat.size1(), j+mat.size2()));
		}
	}

	return mat;
}


/**
 * SU(2) generators, pauli matrices sig_i = 2*S_i
 * @see (Arfken 2013), p. 110
 */
template<class t_mat>
const t_mat& su2_matrix(std::size_t which)
requires is_mat<t_mat> && is_complex<typename t_mat::value_type>
{
	using t_cplx = typename t_mat::value_type;
	constexpr t_cplx c0(0,0);
	constexpr t_cplx c1(1,0);
	constexpr t_cplx cI(0,1);

	static const t_mat mat[] =
	{
		create<t_mat>({{c0, c1}, { c1,  c0}}),	// x
		create<t_mat>({{c0, cI}, {-cI,  c0}}),	// y
		create<t_mat>({{c1, c0}, { c0, -c1}}),	// z
	};

	return mat[which];
}


/**
 * get a vector of pauli matrices
 * @see (Arfken 2013), p. 110
 */
template<class t_vec>
t_vec su2_matrices(bool bIncludeUnit = false)
requires is_basic_vec<t_vec> && is_mat<typename t_vec::value_type>
	&& is_complex<typename t_vec::value_type::value_type>
{
	using t_mat = typename t_vec::value_type;

	t_vec vec;
	vec.reserve(4);

	if(bIncludeUnit)
		vec.emplace_back(unit<t_mat>(2));
	for(std::size_t i=0; i<3; ++i)
		vec.emplace_back(su2_matrix<t_mat>(i));

	return vec;
}


/**
 * project the vector of SU(2) matrices onto a vector
 * proj = <sigma|vec>
 */
template<class t_vec, class t_mat>
t_mat proj_su2(const t_vec& vec, bool is_normalised=1)
requires is_vec<t_vec> && is_mat<t_mat>
{
	typename t_vec::value_type len = 1;
	if(!is_normalised)
		len = norm<t_vec>(vec);

	const auto sigma = su2_matrices<std::vector<t_mat>>(false);
	return inner<std::vector<t_mat>, t_vec>(sigma, vec);
}


/**
 * SU(2) ladders
 * @see https://en.wikipedia.org/wiki/Ladder_operator
 */
template<class t_mat>
const t_mat& su2_ladder(std::size_t which)
requires is_mat<t_mat> && is_complex<typename t_mat::value_type>
{
	using t_cplx = typename t_mat::value_type;
	constexpr t_cplx cI(0,1);
	constexpr t_cplx c05(0.5, 0);

	static const t_mat mat[] =
	{
		c05*su2_matrix<t_mat>(0) + c05*cI*su2_matrix<t_mat>(1),	// up
		c05*su2_matrix<t_mat>(0) - c05*cI*su2_matrix<t_mat>(1),	// down
	};

	return mat[which];
}


/**
 * SU(3) generators, Gell-Mann matrices
 * @see https://de.wikipedia.org/wiki/Gell-Mann-Matrizen
 */
template<class t_mat>
const t_mat& su3_matrix(std::size_t which)
requires is_mat<t_mat> && is_complex<typename t_mat::value_type>
{
	using t_cplx = typename t_mat::value_type;
	using t_real = typename t_cplx::value_type;
	constexpr t_cplx c0(0,0);
	constexpr t_cplx c1(1,0);
	constexpr t_cplx c2(2,0);
	constexpr t_cplx cI(0,1);
	constexpr t_real s3 = std::sqrt(3.);

	static const t_mat mat[] =
	{
		create<t_mat>({{c0,c1,c0}, {c1,c0,c0}, {c0,c0,c0}}),			// 1
		create<t_mat>({{c0,cI,c0}, {-cI,c0,c0}, {c0,c0,c0}}),			// 2
		create<t_mat>({{c1,c0,c0}, {c0,-c1,c0}, {c0,c0,c0}}),			// 3
		create<t_mat>({{c0,c0,c1}, {c0,c0,c0}, {c1,c0,c0}}),			// 4
		create<t_mat>({{c0,c0,cI}, {c0,c0,c0}, {-cI,c0,c0}}),			// 5
		create<t_mat>({{c0,c0,c0}, {c0,c0,c1}, {c0,c1,c0}}),			// 6
		create<t_mat>({{c0,c0,c0}, {c0,c0,cI}, {c0,-cI,c0}}),			// 7
		create<t_mat>({{c1/s3,c0,c0}, {c0,c1/s3,c0}, {c0,c0,-c2/s3*c1}}),	// 8
	};

	return mat[which];
}


/**
 * real crystallographic A matrix
 * @see https://en.wikipedia.org/wiki/Fractional_coordinates
 *
 * derivation from dot products:
 *   (1) <a|b> = ab cos(_cc), (2) <a|c> = ac cos(_bb), (3) <b|c> = bc cos(_aa),
 *   (4) <a|a> = a^2,         (5) <b|b> = b^2,         (6) <c|c> = c^2
 *
 *   (1) => b_0 = b cos(_cc)
 *   (5) => b_1 = b sin(_cc)

 *   (2) => c_0 = c cos(_bb)
 *   (3) => b_0*c_0 + b_1*c_1 = bc cos(_aa)
 *       => bc cos(_cc)*cos(_bb) + b sin(_cc) c_1 = bc cos(_aa)
 *       => c_1 = (c cos(_aa) - c cos(_cc)*cos(_bb)) / sin(_cc)
 *   (6) => c_2^2 = c^2 - c_0^2 - c_1^2
 *       => c_2^2 = c^2 (sin^2(_bb) - (cos(_aa) - cos(_cc)*cos(_bb))^2 / sin^2(_cc))
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat A_matrix(t_real a, t_real b, t_real c, t_real _aa, t_real _bb, t_real _cc)
requires is_mat<t_mat>
{
	const t_real ca = std::cos(_aa);
	const t_real cb = std::cos(_bb);
	const t_real cc = std::cos(_cc);
	const t_real sc = std::sin(_cc);
	const t_real sb = std::sin(_bb);

	return create<t_mat>({
		a,         b*cc,      c*cb,
		t_real{0}, b*sc,      c*(ca - cc*cb)/sc,
		t_real{0}, t_real{0}, c*std::sqrt(sb*sb - std::pow((ca - cc*cb)/sc, t_real{2}))
	});
}


/**
 * reciprocal crystallographic B matrix, B = 2pi * A^(-T)
 * @see https://en.wikipedia.org/wiki/Fractional_coordinates
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat B_matrix(t_real a, t_real b, t_real c, t_real _aa, t_real _bb, t_real _cc)
requires is_mat<t_mat>
{
	const t_real sc = std::sin(_cc);
	const t_real ca = std::cos(_aa);
	const t_real cb = std::cos(_bb);
	const t_real cc = std::cos(_cc);
	const t_real rr = std::sqrt(
		t_real{1} + t_real{2}*ca*cb*cc - (ca*ca + cb*cb + cc*cc));

	return t_real{2}*pi<t_real> * create<t_mat>({
		t_real{1}/a,             t_real{0},                   t_real{0},
		-t_real{1}/a * cc/sc,    t_real{1}/b * t_real{1}/sc,  t_real{0},
		(cc*ca - cb)/(a*sc*rr),  (cb*cc - ca)/(b*sc*rr),      sc/(c*rr)
	});
}


/**
 * UB orientation matrix
 * @see https://dx.doi.org/10.1107/S0021889805004875
 */
template<class t_mat, class t_vec>
t_mat UB_matrix(const t_mat& B,
	const t_vec& vec1_rlu, const t_vec& vec2_rlu, const t_vec& vec3_rlu)
requires is_mat<t_mat> && is_vec<t_vec>
{
	t_vec vec1_lab = B * vec1_rlu;
	t_vec vec2_lab = B * vec2_rlu;
	t_vec vec3_lab = B * vec3_rlu;

	vec1_lab /= norm<t_vec>(vec1_lab);
	vec2_lab /= norm<t_vec>(vec2_lab);
	vec3_lab /= norm<t_vec>(vec3_lab);

	t_mat U_lab = unit<t_mat>(B.size1(), B.size2());
	set_row<t_mat, t_vec>(U_lab, vec1_lab, 0);
	set_row<t_mat, t_vec>(U_lab, vec2_lab, 1);
	set_row<t_mat, t_vec>(U_lab, vec3_lab, 2);

	return U_lab * B;
}


/**
 * get correct distance in unit cell, considering wrapping-around
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
t_real get_dist_uc(const t_mat& matA, const t_vec& vec1, const t_vec& vec2)
requires is_mat<t_mat> && is_vec<t_vec>
{
	t_vec vec1A = matA * vec1;

	// all supercell position to try
	std::vector<t_vec> vecSCs
	{{
		create<t_vec>({0,  0,  0}),

		create<t_vec>({+1,  0,  0}), create<t_vec>({ 0, +1,  0}), create<t_vec>({ 0,  0, +1}),

		create<t_vec>({ 0, +1, +1}), create<t_vec>({ 0, +1, -1}),
		create<t_vec>({+1,  0, +1}), create<t_vec>({+1,  0, -1}),
		create<t_vec>({+1, +1,  0}), create<t_vec>({+1, -1,  0}),

		create<t_vec>({+1, +1, +1}), create<t_vec>({+1, +1, -1}),
		create<t_vec>({+1, -1, +1}), create<t_vec>({+1, -1, -1}),
	}};

	t_real thedist = std::numeric_limits<t_real>::max();

	for(const t_vec& _vecSC : vecSCs)
	{
		for(const t_vec& vecSC : { _vecSC, -_vecSC, t_real{2}*_vecSC, t_real{-2}*_vecSC })
		{
			t_vec vec2A = matA * (vec2 + vecSC);
			t_real dist = norm(vec1A - vec2A);

			thedist = std::min(thedist, dist);
		}
	}

	return thedist;
}


/**
 * conjugate complex vector
 */
template<class t_vec>
t_vec conj(const t_vec& vec)
requires is_basic_vec<t_vec>
{
	const std::size_t N = vec.size();
	t_vec vecConj = zero<t_vec>(N);

	for(std::size_t iComp=0; iComp<N; ++iComp)
	{
		if constexpr(is_complex<typename t_vec::value_type>)
			vecConj[iComp] = std::conj(vec[iComp]);
		else	// simply copy non-complex vector
			vecConj[iComp] = vec[iComp];
	}

	return vecConj;
}


/**
 * conjugate complex matrix
 */
template<class t_mat>
t_mat conj(const t_mat& mat)
requires is_basic_mat<t_mat>
{
	t_mat mat2;
	if constexpr(is_dyn_mat<t_mat>)
		mat2 = t_mat(mat.size2(), mat.size1());

	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			if constexpr(is_complex<typename t_mat::value_type>)
				mat2(i,j) = std::conj(mat(i,j));
			else	// simply transpose non-complex matrix
				mat2(i,j) = mat(i,j);
		}
	}

	return mat2;
}


/**
 * hermitian conjugate complex matrix
 */
template<class t_mat>
t_mat herm(const t_mat& mat)
requires is_basic_mat<t_mat>
{
	t_mat mat2;
	if constexpr(is_dyn_mat<t_mat>)
		mat2 = t_mat(mat.size2(), mat.size1());

	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			if constexpr(is_complex<typename t_mat::value_type>)
				mat2(j,i) = std::conj(mat(i,j));
			else	// simply transpose non-complex matrix
				mat2(j,i) = mat(i,j);
		}
	}

	return mat2;
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
/**
 * wrap atom positions back to unit cell
 */
template<class t_vec, class t_real = typename t_vec::value_type>
t_vec keep_atom_in_uc(const t_vec& _atom, std::size_t dim = 3)
requires is_vec<t_vec>
{
	t_vec newatom = _atom;

	for(std::size_t i=0; i</*newatom.size()*/ dim; ++i)
	{
		newatom[i] = std::fmod(newatom[i], t_real{1});
		if(newatom[i] <= t_real{-0.5})
			newatom[i] += std::abs(std::floor(newatom[i]));
		if(newatom[i] > t_real{0.5})
			newatom[i] -= std::abs(std::ceil(newatom[i]));
	}

	return newatom;
}



/**
 * wrap collection of atom positions back to unit cell
 */
template<class t_vec, class t_real = typename t_vec::value_type,
	template<class...> class t_cont = std::vector>
t_cont<t_vec> keep_atoms_in_uc(const t_cont<t_vec>& _atoms, std::size_t dim = 3)
requires is_vec<t_vec>
{
	t_cont<t_vec> newatoms;
	newatoms.reserve(_atoms.size());

	for(const auto& _atom : _atoms)
		newatoms.emplace_back(keep_atom_in_uc<t_vec, t_real>(_atom, dim));

	return newatoms;
}



/**
 * get the supercell vector and the unitcell atom index from a supercell atom position
 */
template<class t_vec, class t_real = typename t_vec::value_type,
	template<class...> class t_cont = std::vector>
std::tuple<bool, std::size_t, t_vec>
get_supercell(const t_vec& sc_pos, const t_cont<t_vec>& _uc_sites,
	std::size_t dim = 3,
	t_real eps = std::numeric_limits<t_real>::eps())
requires is_vec<t_vec>
{
	t_cont<t_vec> uc_sites = keep_atoms_in_uc<t_vec, t_real, t_cont>(_uc_sites, dim);

	// unit cell position corresponding to supercell position
	t_vec uc_pos = keep_atom_in_uc<t_vec, t_real>(sc_pos, dim);

	// get corresponding unit cell index
	bool found = false;
	std::size_t uc_site_idx = 0;
	t_vec sc_vec = sc_pos - uc_pos;

	if(auto iter = std::find_if(uc_sites.begin(), uc_sites.end(),
		[&uc_pos, &eps](const t_vec& vec) -> bool
	{
		return equals<t_vec>(uc_pos, vec, eps);
	}); iter != uc_sites.end())
	{
		found = true;
		uc_site_idx = iter - uc_sites.begin();
	}

	return std::make_tuple(found, uc_site_idx, sc_vec);
}



/**
 * create positions using the given symmetry operations
 */
template<class t_vec, class t_mat, class t_real = typename t_vec::value_type,
	template<class...> class t_cont = std::vector>
t_cont<t_vec> apply_ops_hom(const t_vec& _atom, const t_cont<t_mat>& ops,
	t_real eps=std::numeric_limits<t_real>::epsilon(),
	bool keepInUnitCell = true, bool ignore_occupied = false, bool ret_hom = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	// in homogeneous coordinates
	t_vec atom = _atom;
	if(atom.size() == 3)
		atom = create<t_vec>({atom[0], atom[1], atom[2], 1});

	t_cont<t_vec> newatoms;
	newatoms.reserve(ops.size());

	for(const auto& op : ops)
	{
		auto newatom = op*atom;

		// return a homogeneous 4-vector or a 3-vector
		if(!ret_hom)
			newatom.resize(3);

		if(keepInUnitCell)
			newatom = keep_atom_in_uc<t_vec>(newatom, 3);

		// position already occupied?
		if(ignore_occupied ||
			std::find_if(newatoms.begin(), newatoms.end(),
				[&newatom, eps](const t_vec& vec) -> bool
		{
			return tl2::equals<t_vec>(vec, newatom, eps);
		}) == newatoms.end())
		{
			newatoms.emplace_back(std::move(newatom));
		}
	}

	return newatoms;
}
// ----------------------------------------------------------------------------

}



// ----------------------------------------------------------------------------
// lapack wrappers
// ----------------------------------------------------------------------------
#ifdef __TLIBS2_USE_LAPACK__

namespace tl2_la {

/**
 * LU decomposition of a matrix, mat = P * L * U, returning raw results
 * @returns [ok, LU, perm]
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dgetrf.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zgetrf.html
 */
template<class t_mat, template<class...> class t_vec = std::vector>
std::tuple<bool, t_vec<typename t_mat::value_type>, t_vec<lapack_int>> _lu_raw(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;
	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();
	const std::size_t minor = std::min(rows, cols);


	t_vec<t_scalar> outmat(rows*cols);
	t_vec<lapack_int> outpivots(minor);

	for(std::size_t i=0; i<rows; ++i)
		for(std::size_t j=0; j<cols; ++j)
			outmat[i*cols + j] = mat(i, j);

	int err = -1;
	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cgetrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outpivots.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outpivots.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_sgetrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outpivots.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outpivots.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}

	return std::make_tuple(err == 0, outmat, outpivots);
}


/**
 * LU decomposition of a matrix, mat = P * L * U
 * @returns [ok, P, L, U]
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dgetrf.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zgetrf.html
 */
template<class t_mat, template<class...> class t_vec = std::vector>
std::tuple<bool, t_mat, t_mat, t_mat> lu(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;

	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();


	auto [ ok, lumat, pivots ] = _lu_raw<t_mat, t_vec>(mat);

	t_mat P = tl2::unit<t_mat>(rows, cols);
	t_mat L = tl2::unit<t_mat>(rows, cols);
	t_mat U = tl2::unit<t_mat>(rows, cols);

	// L and U
	for(std::size_t i=0; i<rows; ++i)
	{
		for(std::size_t j=0; j<cols; ++j)
		{
			if(j>=i)
				U(i, j) = lumat[i*cols + j];
			else
				L(i, j) = lumat[i*cols + j];
		}
	}

	// permutation matrix P
	for(std::size_t i=0; i<pivots.size(); ++i)
		P = tl2::prod<t_mat>(P, tl2::perm<t_mat>(rows, cols, i, pivots[i]-1));

	return std::make_tuple(ok, P, L, U);
}


/**
 * inverted matrix
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dgetri.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zgetri.html
 */
template<class t_mat>
std::tuple<t_mat, bool> inv(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	// fail if matrix is not square
	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size1() == mat.size2()));
	else
		static_assert(t_mat::size1() == t_mat::size2());

	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();

	t_mat I = tl2::unit<t_mat>(rows, cols);


	// lu factorisation
	auto [ ok, lumat, pivots ] = _lu_raw<t_mat, std::vector>(mat);
	if(!ok)
		return std::make_tuple(I, false);


	// inversion
	int err = -1;
	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cgetri(LAPACK_ROW_MAJOR,
				rows, lumat.data(), rows,
				pivots.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zgetri(LAPACK_ROW_MAJOR,
				rows, lumat.data(), rows,
				pivots.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_sgetri(LAPACK_ROW_MAJOR,
				rows, lumat.data(), rows,
				pivots.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dgetri(LAPACK_ROW_MAJOR,
				rows, lumat.data(), rows,
				pivots.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}


	for(std::size_t i=0; i<rows; ++i)
		for(std::size_t j=0; j<cols; ++j)
			I(i, j) = lumat[i*cols + j];

	return std::make_tuple(I, err == 0);
}


/**
 * QR decomposition of a matrix, mat = QR
 * @returns [ok, Q, R]
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dgeqrf.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zgeqrf.html
 */
template<class t_mat, class t_vec = std::vector<typename t_mat::value_type>>
std::tuple<bool, t_mat, t_mat> qr(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;
	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();
	const std::size_t minor = std::min(rows, cols);

	const t_mat I = tl2::unit<t_mat>(minor);
	t_mat Q = I, R = mat;


	t_vec outmat(rows*cols), outvec(minor);

	for(std::size_t i=0; i<rows; ++i)
		for(std::size_t j=0; j<cols; ++j)
			outmat[i*cols + j] = mat(i, j);

	int err = -1;
	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cgeqrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outvec.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zgeqrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outvec.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_sgeqrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outvec.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outvec.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}

	for(std::size_t i=0; i<rows; ++i)
		for(std::size_t j=0; j<cols; ++j)
			R(i, j) = (j>=i ? outmat[i*cols + j] : t_real{0});


	t_vec v = tl2::zero<t_vec>(minor);

	for(std::size_t k=1; k<=minor; ++k)
	{
		for(std::size_t i=1; i<=k-1; ++i)
			v[i-1] = t_real{0};
		v[k-1] = t_real{1};

		for(std::size_t i=k+1; i<=minor; ++i)
			v[i-1] = outmat[(i-1)*cols + (k-1)];

		Q = Q * (I - outvec[k-1]*tl2::outer<t_mat, t_vec>(v, v));
	}

	return std::make_tuple(err == 0, Q, R);
}


/**
 * cholesky decomposition of a hermitian, positive-definite matrix, mat = C^H C
 * @returns [ok, C]
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dpotrf.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zpotrf.html
 * @see http://www.netlib.org/utk/papers/factor/node9.html
 */
template<class t_mat, class t_vec = std::vector<typename t_mat::value_type>>
std::tuple<bool, t_mat> chol(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;
	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size1() == mat.size2()));
	else
		static_assert(t_mat::size1() == t_mat::size2());

	const std::size_t N = mat.size1();
	int err = -1;

	std::vector<t_scalar> outmat(N*N);

	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<N; ++j)
			outmat[i*N + j] = (j >= i ? mat(i, j) : 0.);

	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cpotrf(LAPACK_ROW_MAJOR,
				'U', N, outmat.data(), N);
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zpotrf(LAPACK_ROW_MAJOR,
				'U', N, outmat.data(), N);
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_spotrf(LAPACK_ROW_MAJOR,
				'U', N, outmat.data(), N);
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dpotrf(LAPACK_ROW_MAJOR,
				'U', N, outmat.data(), N);
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}

	t_mat C = tl2::create<t_mat>(N, N);

	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<N; ++j)
			C(i, j) = outmat[i*N + j];

	//std::cerr << "error value: " << err << std::endl;
	return std::make_tuple(err == 0, C);
}


/**
 * cholesky decomposition of a hermitian matrix, mat = C D C^H
 * @returns [ok, C, D]
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zhptrf.html
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dsptrf.html
 */
template<class t_mat, class t_vec = std::vector<typename t_mat::value_type>>
std::tuple<bool, t_mat, t_mat> chol2(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;
	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size1() == mat.size2()));
	else
		static_assert(t_mat::size1() == t_mat::size2());

	const std::size_t N = mat.size1();
	int err = -1;

	std::vector<t_scalar> outmat(N*(N+1)/2);
	std::vector<int> outvec(N);

	std::size_t lin_idx = 0;
	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<=i; ++j)
			outmat[lin_idx++] = mat(i, j);

	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_chptrf(LAPACK_ROW_MAJOR,
				'L', N, outmat.data(), outvec.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zhptrf(LAPACK_ROW_MAJOR,
				'L', N, outmat.data(), outvec.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_ssptrf(LAPACK_ROW_MAJOR,
				'L', N, outmat.data(), outvec.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dsptrf(LAPACK_ROW_MAJOR,
				'L', N, outmat.data(), outvec.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}

	t_mat C = tl2::unit<t_mat>(N, N);
	t_mat D = tl2::unit<t_mat>(N);

	lin_idx = 0;
	for(std::size_t i=0; i<N; ++i)
	{
		for(std::size_t j=0; j<=i; ++j)
		{
			if(i == j)
				D(i, j) = outmat[lin_idx++];
			else
				C(i, j) = outmat[lin_idx++];
		}
	}

	// TODO: only correct for unit permutation so far...
	for(int& i : outvec)
		i = std::abs(i)-1;
	t_mat perm = tl2::perm<t_mat>(outvec);
	C = perm * C;

	//std::cerr << "error value: " << err << std::endl;
	return std::make_tuple(err == 0, C, D);
}


/**
 * eigenvectors and -values of a complex matrix
 * @returns [ok, evals, evecs]
 */
template<class t_mat_cplx, class t_vec_cplx, class t_cplx = typename t_mat_cplx::value_type,
    class t_real = typename t_cplx::value_type>
std::tuple<bool, std::vector<t_cplx>, std::vector<t_vec_cplx>>
eigenvec(const t_mat_cplx& mat,
	bool only_evals=false, bool is_hermitian=false, bool normalise=false,
	t_real mineval=-1, t_real maxeval=-2, t_real eps=-1)
	requires tl2::is_mat<t_mat_cplx> && tl2::is_vec<t_vec_cplx> && tl2::is_complex<t_cplx>
{
	bool only_selected_evals = (mineval <= maxeval);
	bool use_selective_func = only_selected_evals;
	//use_selective_func = true;

	std::vector<t_cplx> evals;
	std::vector<t_vec_cplx> evecs;

	if(mat.size1() != mat.size2() || mat.size1() == 0)
		return std::make_tuple(0, evals, evecs);

	const std::size_t N = mat.size1();
	evals.resize(N);

	if(!only_evals)
		evecs.resize(N, tl2::zero<t_vec_cplx>(N));

	std::vector<t_cplx> inmat(N*N, t_cplx{0,0}),
		outevecs(only_evals ? 0 : N*N, t_cplx{0,0});

	for(std::size_t i=0; i<N; ++i)
	{
		for(std::size_t j=0; j<N; ++j)
		{
			if(is_hermitian)
				inmat[i*N + j] = (j>=i ? mat(j,i) : t_real{0});
			else
				inmat[i*N + j] = mat(j,i);
		}
	}


	int err = -1;

	if(is_hermitian)
	{
		// evals of hermitian matrix are purely real
		std::vector<t_real> outevals_real(N, t_real{0});

		// all eigenvalues
		if(!use_selective_func)
		{
			if constexpr(std::is_same_v<t_real, float>)
			{
				err = LAPACKE_cheev(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					'L', N, inmat.data(), N,
					outevals_real.data());
			}
			else if constexpr(std::is_same_v<t_real, double>)
			{
				err = LAPACKE_zheev(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					'L', N, inmat.data(), N,
					outevals_real.data());
			}
			else
			{
				static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
				//throw std::domain_error("Invalid element type.");
			}
		}

		// only selected eigenvalues
		else
		{
			int minidx = 1, maxidx = N;
			int iNumFound = 0;

			std::unique_ptr<int, std::default_delete<int[]>>
			uptrIdxArr(new int[2*N]);

			// use maximum precision if none given
			if(eps < t_real{0})
			{
				if constexpr(std::is_same_v<t_real, float>)
				{
					eps = LAPACKE_slamch('S');
				}
				else if constexpr(std::is_same_v<t_real, double>)
				{
					eps = LAPACKE_dlamch('S');
				}
				else
				{
					static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
					//throw std::domain_error("Invalid element type.");
				}
			}

			if constexpr(std::is_same_v<t_real, float>)
			{
				err = LAPACKE_cheevr(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					only_selected_evals ? 'V' : 'A',
					'L', N, inmat.data(), N,
					mineval, maxeval, minidx, maxidx,
					eps, &iNumFound,
					outevals_real.data(), outevecs.data(),
					N, uptrIdxArr.get());
			}
			else if constexpr(std::is_same_v<t_real, double>)
			{
				err = LAPACKE_zheevr(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					only_selected_evals ? 'V' : 'A',
					'L', N, inmat.data(), N,
					mineval, maxeval, minidx, maxidx,
					eps, &iNumFound,
					outevals_real.data(), outevecs.data(),
					N, uptrIdxArr.get());
			}
			else
			{
				static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
				//throw std::domain_error("Invalid element type.");
			}

			// resize to actual number of eigenvalues and -vectors
			if(std::size_t(iNumFound) != N)
			{
				evals.resize(iNumFound, t_real{0});
				evecs.resize(iNumFound, tl2::zero<t_vec_cplx>(N));
			}
		}

		// copy to complex output vector
		for(std::size_t i=0; i<evals.size(); ++i)
			evals[i] = outevals_real[i];
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cgeev(LAPACK_COL_MAJOR,
				'N', only_evals ? 'N' : 'V', N,
				inmat.data(), N, evals.data(), nullptr, N,
				only_evals ? nullptr : outevecs.data(), N);
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zgeev(LAPACK_COL_MAJOR,
				'N', only_evals ? 'N' : 'V', N,
				inmat.data(), N, evals.data(), nullptr, N,
				only_evals ? nullptr : outevecs.data(), N);
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid real type.");
		}
	}


	if(!only_evals)
	{
		for(std::size_t i=0; i<evecs.size(); ++i)
		{
			// hermitian algo overwrites original matrix!
			for(std::size_t j=0; j<N; ++j)
				evecs[i][j] = (is_hermitian && !use_selective_func)
					? inmat[i*N + j] : outevecs[i*N + j];

			if(normalise && (err == 0))
			{
				t_cplx n = tl2::norm(evecs[i]);
				if(!tl2::equals<t_cplx>(n, t_cplx{0,0}))
					evecs[i] /= n;
			}
		}
	}

	return std::make_tuple(err == 0, evals, evecs);
}


/**
 * eigenvectors and -values of a real matrix
 * @returns [ok, evals_re, evals_im, evecs_re, evecs_im]
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
std::tuple<bool, std::vector<t_real>, std::vector<t_real>, std::vector<t_vec>, std::vector<t_vec>>
eigenvec(const t_mat& mat, bool only_evals=false, bool is_symmetric=false, bool normalise=false,
	t_real mineval=-1, t_real maxeval=-2, t_real eps=-1)
	requires (tl2::is_mat<t_mat> && tl2::is_vec<t_vec> && !tl2::is_complex<t_real>)
{
	bool only_selected_evals = (mineval <= maxeval);
	bool use_selective_func = only_selected_evals;
	//use_selective_func = true;

	std::vector<t_real> evals_re, evals_im;
	std::vector<t_vec> evecs_re, evecs_im;

	if(mat.size1() != mat.size2() || mat.size1() == 0)
		return std::make_tuple(false, evals_re, evals_im, evecs_re, evecs_im);

	const std::size_t N = mat.size1();
	evals_re.resize(N, t_real{0});
	evals_im.resize(N, t_real{0});

	if(!only_evals)
	{
		evecs_re.resize(N, tl2::zero<t_vec>(N));
		evecs_im.resize(N, tl2::zero<t_vec>(N));
	}

	std::vector<t_real> inmat(N*N, t_real{0}),
		outevecs(only_evals ? 0 : N*N, t_real{0});

	for(std::size_t i=0; i<N; ++i)
	{
		for(std::size_t j=0; j<N; ++j)
		{
			if(is_symmetric)
				inmat[i*N + j] = (j>=i ? mat(j,i) : t_real{0});
			else
				inmat[i*N + j] = mat(j,i);
		}
	}

	int err = -1;

	if(is_symmetric)
	{
		// all eigenvalues
		if(!use_selective_func)
		{
			// evals of symmetric matrix are purely real
			if constexpr(std::is_same_v<t_real, float>)
			{
				err = LAPACKE_ssyev(LAPACK_COL_MAJOR,
					(only_evals ? 'N' : 'V'),
					'L', N, inmat.data(), N,
					evals_re.data());
			}
			else if constexpr(std::is_same_v<t_real, double>)
			{
				err = LAPACKE_dsyev(LAPACK_COL_MAJOR,
					(only_evals ? 'N' : 'V'),
					'L', N, inmat.data(), N,
					evals_re.data());
			}
			else
			{
				static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
				//throw std::domain_error("Invalid element type.");
			}
		}

		// only selected eigenvalues
		else
		{
			int minidx = 1, maxidx = N;
			int iNumFound = 0;

			std::unique_ptr<int, std::default_delete<int[]>>
			uptrIdxArr(new int[2*N]);

			// use maximum precision if none given
			if(eps < t_real{0})
			{
				if constexpr(std::is_same_v<t_real, float>)
				{
					eps = LAPACKE_slamch('S');
				}
				else if constexpr(std::is_same_v<t_real, double>)
				{
					eps = LAPACKE_dlamch('S');
				}
				else
				{
					static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
					//throw std::domain_error("Invalid element type.");
				}
			}

			if constexpr(std::is_same_v<t_real, float>)
			{
				err = LAPACKE_ssyevr(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					only_selected_evals ? 'V' : 'A',
					'L', N, inmat.data(), N,
					mineval, maxeval, minidx, maxidx,
					eps, &iNumFound, evals_re.data(),
					outevecs.data(), N, uptrIdxArr.get());
			}
			else if constexpr(std::is_same_v<t_real, double>)
			{
				err = LAPACKE_dsyevr(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					only_selected_evals ? 'V' : 'A',
					'L', N, inmat.data(), N,
					mineval, maxeval, minidx, maxidx,
					eps, &iNumFound, evals_re.data(),
					outevecs.data(), N, uptrIdxArr.get());
			}
			else
			{
				static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
				//throw std::domain_error("Invalid element type.");
			}

			// resize to actual number of eigenvalues and -vectors
			if(std::size_t(iNumFound) != N)
			{
				evals_re.resize(iNumFound, t_real{0});
				evals_im.resize(iNumFound, t_real{0});
				evecs_re.resize(iNumFound, tl2::zero<t_vec>(N));
				evecs_im.resize(iNumFound, tl2::zero<t_vec>(N));
			}
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_sgeev(LAPACK_COL_MAJOR,
				'N', only_evals ? 'N' : 'V', N,
				inmat.data(), N,
				evals_re.data(), evals_im.data(),
				nullptr, N,
				only_evals ? nullptr : outevecs.data(), N);
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dgeev(LAPACK_COL_MAJOR,
				'N', only_evals ? 'N' : 'V', N,
				inmat.data(), N,
				evals_re.data(), evals_im.data(),
				nullptr, N,
				only_evals ? nullptr : outevecs.data(), N);
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}


	// evecs
	if(!only_evals)
	{
		if((is_symmetric && !use_selective_func))
		{
			for(std::size_t i=0; i<evals_re.size(); ++i)
			{
				// symmetric algo overwrites original matrix!
				for(std::size_t j=0; j<N; ++j)
					evecs_re[i][j] = inmat[i*N + j];
			}
		}
		else
		{
			for(std::size_t i=0; i<evals_re.size(); ++i)
			{
				if(tl2::equals<t_real>(evals_im[i], 0))
				{
					for(std::size_t j=0; j<N; ++j)
					{
						evecs_re[i][j] = outevecs[i*N + j];
						evecs_im[i][j] = 0;
					}
				}
				else
				{
					for(std::size_t j=0; j<N; ++j)
					{
						evecs_re[i][j] = outevecs[i*N + j];
						evecs_im[i][j] = outevecs[(i+1)*N + j];	// imag part of evec follows next in array

						evecs_re[i+1][j] = evecs_re[i][j];		// next evec is the conjugated one
						evecs_im[i+1][j] = -evecs_im[i][j];
					}
					++i;	// already used two values in array
				}
			}
		}

		if(normalise && (err == 0))
		{
			for(std::size_t i=0; i<evecs_re.size(); ++i)
			{
				t_real sum{0};
				for(std::size_t j=0; j<N; ++j)
					sum += std::norm(
						std::complex(evecs_re[i][j],
							evecs_im[i][j]));
				sum = std::sqrt(sum);

				if(!tl2::equals<t_real>(sum, 0))
				{
					evecs_re[i] /= sum;
					evecs_im[i] /= sum;
				}
			}
		}
	}

	return std::make_tuple(err == 0, evals_re, evals_im, evecs_re, evecs_im);
}


/**
 * singular values of a real or complex matrix mat = U * diag{vals} * V^h
 * @returns [ ok, U, Vh, vals ]
 */
template<class t_mat, class t_scalar = typename t_mat::value_type, class t_real = tl2::underlying_value_type<t_scalar>>
std::tuple<bool, t_mat, t_mat, std::vector<t_real>>
singval(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();

	const auto [Nmin, Nmax] = std::minmax(rows, cols);
	std::vector<t_scalar> inmat(rows*cols), outU(rows*rows), outVh(cols*cols);
	std::vector<t_real> vals(Nmin);
	std::vector<t_real> _tmp(Nmax * Nmax * 2);

	for(std::size_t i=0; i<rows; ++i)
		for(std::size_t j=0; j<cols; ++j)
			inmat[i*cols + j] = mat(i,j);

	int err = -1;
	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cgesvd(LAPACK_ROW_MAJOR,
				'A', 'A', rows, cols,
				inmat.data(), cols,
				vals.data(), outU.data(),
				rows, outVh.data(),
				cols, _tmp.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zgesvd(LAPACK_ROW_MAJOR,
				'A', 'A', rows, cols,
				inmat.data(), cols,
				vals.data(), outU.data(),
				rows, outVh.data(),
				cols, _tmp.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_sgesvd(LAPACK_ROW_MAJOR,
				'A', 'A', rows, cols, inmat.data(), cols,
				vals.data(), outU.data(),
				rows, outVh.data(), cols, _tmp.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,
				'A', 'A', rows, cols, inmat.data(), cols,
				vals.data(), outU.data(),
				rows, outVh.data(), cols, _tmp.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			//throw std::domain_error("Invalid element type.");
		}
	}


	t_mat U = tl2::unit<t_mat>(rows);
	t_mat Vh = tl2::unit<t_mat>(cols);

	for(std::size_t i=0; i<Nmax; ++i)
	{
		for(std::size_t j=0; j<Nmax; ++j)
		{
			if(i<U.size1() && j<U.size2()) U(i,j) = outU[i*cols + j];
			if(i<Vh.size1() && j<Vh.size2()) Vh(i,j) = outVh[i*cols + j];
		}
	}

	return std::make_tuple(err==0, U, Vh, vals);
}


/**
 * pseudoinverse M+ of a matrix
 * @see https://de.wikipedia.org/wiki/Pseudoinverse#Berechnung
 * @see (Arens 2015), pp. 788-792
 *
 * M  = U D (V*)^h
 * M+ = V D+ (U*)^h
 */
template<class t_mat>
std::tuple<t_mat, bool> pseudoinv(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	auto [ ok, U, Vh, vals ] = singval<t_mat>(mat);

	auto V = tl2::herm(Vh);
	auto Uh = tl2::herm(U);

	for(t_real& d : vals)
	{
		if(!tl2::equals<t_real>(d, t_real(0)))
			d = t_real(1)/d;
	}

	auto diag = tl2::diag<t_mat>(vals);
	return std::make_tuple(tl2::prod<t_mat>(V, tl2::prod(diag, Uh)), ok);
}



// ----------------------------------------------------------------------------
// equation solvers
// ----------------------------------------------------------------------------

/**
 * system of ODEs with constant coefficients C
 * @see (Arens 2015), pp. 1049-1051
 *
 * f'(x) = C f(x) and f(x0) = f0
 * => f(x) = f0 * exp(C(x-x0)) = sum_i norm_i * evec_i * exp(eval_i * (x-x0))
 *    norm = (evec_i)^(-1) * f0
 */
template<class t_mat, class t_vec, class t_val = typename t_vec::value_type>
std::tuple<bool, t_vec>
odesys_const(const t_mat& C, const t_val& x, const t_val& x0, const t_vec& f0)
requires (tl2::is_mat<t_mat> && tl2::is_vec<t_vec>)
{
	const std::size_t rows = C.size1();
	const std::size_t cols = C.size2();
	if(rows != cols)
		return std::make_tuple(false, t_vec{});

	const auto [ok, evals, evecs] = eigenvec<t_mat, t_vec, t_val>(C);
	if(!ok)
		return std::make_tuple(false, t_vec{});

	//for(std::size_t i=0; i<evals.size(); ++i)
	//	std::cout << "eval: " << evals[i] << ", evec: " << evecs[i] << std::endl;

	const t_mat basis = tl2::create<t_mat, t_vec, std::vector>(evecs, false);
	const auto [basis_inv, invok] = tl2::inv<t_mat>(basis);
	if(!invok)
		return std::make_tuple(false, t_vec{});
	const t_vec norm = basis_inv * f0;

	t_vec f = tl2::zero<t_vec>(cols);
	for(std::size_t i=0; i<cols; ++i)
		f += norm[i] * evecs[i] * std::exp(evals[i]*(x-x0));

	return std::make_tuple(true, f);
}


/**
 * system of difference equations with constant coefficients C
 */
template<class t_mat, class t_vec, class t_val = typename t_vec::value_type>
std::tuple<bool, t_vec>
diffsys_const(const t_mat& C, const t_val& n, const t_val& n0, const t_vec& f0)
requires (tl2::is_mat<t_mat> && tl2::is_vec<t_vec>)
{
	const std::size_t rows = C.size1();
	const std::size_t cols = C.size2();
	if(rows != cols)
		return std::make_tuple(false, t_vec{});

	const auto [ok, evals, evecs] = eigenvec<t_mat, t_vec, t_val>(C);
	if(!ok)
		return std::make_tuple(false, t_vec{});

	//for(std::size_t i=0; i<evals.size(); ++i)
	//	std::cout << "eval: " << evals[i] << ", evec: " << evecs[i] << std::endl;

	const t_mat basis = tl2::create<t_mat, t_vec, std::vector>(evecs, false);
	const auto [basis_inv, invok] = tl2::inv<t_mat>(basis);
	if(!invok)
		return std::make_tuple(false, t_vec{});
	const t_vec norm = basis_inv * f0;

	t_vec f = tl2::zero<t_vec>(cols);
	for(std::size_t i=0; i<cols; ++i)
		f += norm[i] * evecs[i] * std::pow(evals[i], n-n0);

	return std::make_tuple(true, f);
}
// ----------------------------------------------------------------------------

}

#endif	// __TLIBS2_USE_LAPACK__


namespace tl2 {

/**
 * QR decomposition of a matrix
 * @returns [ok, Q, R]
 * @see (Scarpino 2011), pp. 269-272
 */
template<class t_mat, class t_vec>
std::tuple<bool, t_mat, t_mat> qr(const t_mat& mat)
requires is_mat<t_mat> && is_vec<t_vec>
{
#ifdef __TLIBS2_USE_LAPACK__
	return tl2_la::qr<t_mat, t_vec>(mat);

#else
	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();
	const std::size_t N = std::min(cols, rows);

#if __TLIBS2_QR_METHOD == 0
	t_mat R = mat;
	t_mat Q = unit<t_mat>(N, N);

	for(std::size_t icol=0; icol<N-1; ++icol)
	{
		t_vec vecCol = col<t_mat, t_vec>(R, icol);
		t_mat matMirror = ortho_mirror_zero_op<t_mat, t_vec>(vecCol, icol);

		Q = prod(Q, matMirror, false);
		R = prod(matMirror, R);
	}

#elif __TLIBS2_QR_METHOD == 1
	std::vector<t_vec> sysM;
	sysM.reserve(mat.size2());
	for(std::size_t i=0; i<mat.size2(); ++i)
		sysM.push_back(col<t_mat, t_vec>(mat, i));

	std::vector<t_vec> Qsys = orthonorm_sys<t_vec, std::vector, std::vector>(sysM);
	t_mat Q = create<t_mat>(Qsys, false);
	t_mat R = prod(trans(Q), mat);
#endif

	return std::make_tuple(true, Q, R);
#endif
}


/**
 * inverted matrix
 * @see https://en.wikipedia.org/wiki/Invertible_matrix#In_relation_to_its_adjugate
 * @see https://en.wikipedia.org/wiki/Adjugate_matrix
 */
template<class t_mat>
std::tuple<t_mat, bool> inv(const t_mat& mat)
requires is_mat<t_mat>
{
	// fails if matrix is not square, TODO: fix
	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size1() == mat.size2()));
	else
		static_assert(t_mat::size1() == t_mat::size2());

#ifdef __TLIBS2_USE_LAPACK__
	return tl2_la::inv<t_mat>(mat);
#else
	using T = typename t_mat::value_type;
	using t_vec = std::vector<T>;
	const std::size_t N = mat.size1();

	const auto& matFlat = matvec_adapter<t_mat>{mat};
	const T fullDet = flat_det<t_vec>(matFlat, N);

	// fail if determinant is zero
	if(equals<T>(fullDet, 0))
	{
		//std::cerr << "det == 0" << std::endl;
		return std::make_tuple(t_mat(), false);
	}

	t_mat matInv;
	if constexpr(is_dyn_mat<t_mat>)
		matInv = t_mat(N, N);

	for(std::size_t i=0; i<N; ++i)
	{
		for(std::size_t j=0; j<N; ++j)
		{
			const T sgn = ((i+j) % 2) == 0 ? T(1) : T(-1);
			const t_vec subMat = flat_submat<t_vec>(matFlat, N, N, i, j);
			matInv(j,i) = sgn * flat_det<t_vec>(subMat, N-1);
		}
	}

	matInv = matInv / fullDet;
	return std::make_tuple(matInv, true);
#endif
}


/**
 * determinant
 *
 * M v = l v
 * M = V^-1 L V
 * det(M) = det(V^-1 L V) = det(V V^-1 L) = det(L) = l_1*...*l_n
 */
template<class t_mat>
typename t_mat::value_type det(const t_mat& mat)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;
	if(mat.size1() != mat.size2())
		return 0;

#ifdef __TLIBS2_USE_LAPACK__

	using t_vec = vec<T>;

	if constexpr(tl2::is_complex<typename t_mat::value_type>)
	{
		const auto [ok, evals, evecs] =
			tl2_la::eigenvec<t_mat, t_vec, T>(mat, true, false, false);

		T detval{1, 0};
		for(const auto& eval : evals)
			detval *= eval;

		return detval;
	}
	else
	{
		const auto [ok, evalsRe, evalsIm, evecsRe, evecsIm] =
			tl2_la::eigenvec<t_mat, t_vec, T>(mat, true, false, false);

		std::complex<T> detval{1, 0};
		for(std::size_t i=0; i<evalsRe.size(); ++i)
			detval *= std::complex<T>{evalsRe[i], evalsIm[i]};

		return detval.real();
	}

#else

	const auto& matFlat = matvec_adapter<t_mat>{mat};
	return flat_det<std::vector<T>>(matFlat, mat.size1());

#endif
}


/**
 * matrix power
 */
template<class t_mat>
std::tuple<t_mat, bool> pow(const t_mat& mat, int ipow)
requires is_mat<t_mat>
{
	t_mat themat;
	if(mat.size1() != mat.size2())
		return std::make_tuple(themat, false);

	bool ok = true;
	int ipow_pos = ipow<0 ? -ipow : ipow;

	themat = unit<t_mat>(mat.size1());
	for(int i=0; i<ipow_pos; ++i)
		themat = themat*mat;

	if(ipow < 0)
		std::tie(themat, ok) = tl2::inv<t_mat>(themat);

	return std::make_tuple(themat, ok);
}


/**
 * least-squares regression, solves for parameters v
 * @see https://en.wikipedia.org/wiki/Least_squares
 * @see (Arens 2015), p. 793
 *
 * exact equation:
 * 	X v = y
 *
 * approx. equation (normal equation):
 * 	X^t X v = X^t y
 * 	v = inv(X^t X) X^t y
 *
 * 	(QR)^t QR v = X^t y
 * 	R^tQ^t QR v = X^t y
 * 	R^tR v = X^t y
 * 	v = inv(R^tR) X^t y
 */
template<class t_vec, class t_mat = mat<typename t_vec::value_type, std::vector>>
std::tuple<t_vec, bool> leastsq(const t_vec& x, const t_vec& y, std::size_t order,
	bool use_qr=true, bool use_pseudoinv=false)
requires is_vec<t_vec> && is_dyn_mat<t_mat>
{
	// check array sizes, TODO: fix
	if constexpr(tl2::is_dyn_vec<t_vec>)
		assert((x.size() == y.size()));


	using namespace tl2_ops;
	using T = typename t_vec::value_type;

	const std::size_t N = x.size();
	t_mat X{N, order+1};

	for(std::size_t j=0; j<=order; ++j)
		for(std::size_t i=0; i<N; ++i)
			X(i,j) = std::pow(x[i], static_cast<T>(j));


	t_mat XtX;

	if(use_qr)
	{
		auto [ok_qr, Q, R] = qr<t_mat, t_vec>(X);
		if(!ok_qr)
				return std::make_tuple(t_vec{}, false);

		XtX = trans<t_mat>(R) * R;
	}
	else
	{
		XtX = trans<t_mat>(X) * X;
	}

	t_vec Xty = trans<t_mat>(X) * y;
	t_mat Y;
	bool ok = false;

	if(use_pseudoinv)
	{
#ifdef __TLIBS2_USE_LAPACK__
		std::tie(Y, ok) = tl2_la::pseudoinv<t_mat>(XtX);
#else
		#pragma message("leastsq: Pseudo-inverse is not available, using standard inverse instead.")
		std::tie(Y, ok) = inv<t_mat>(XtX);
#endif
	}
	else
	{
		std::tie(Y, ok) = inv<t_mat>(XtX);
	}

	t_vec v = Y * Xty;

	return std::make_tuple(v, ok);
}


/**
 * signed angle between two vectors
 */
template<typename t_vec>
typename t_vec::value_type angle(const t_vec& vec0,
	const t_vec& vec1, const t_vec* pvec_norm = nullptr)
requires is_vec<t_vec>
{
	using namespace tl2_ops;
	using t_real = typename t_vec::value_type;

	if(vec0.size() != vec1.size())
		throw std::runtime_error("angle: vector sizes do not match.");

	if(vec0.size() == 2)
	{
		// signed angles wrt basis
		t_real angle0 = std::atan2(vec0[1], vec0[0]);
		t_real angle1 = std::atan2(vec1[1], vec1[0]);

		return angle1 - angle0;
	}
	if(vec0.size() == 3)
	{
		// cross product gives sine
		t_vec veccross = cross<t_vec>({vec0, vec1});
		t_real dS = norm(veccross);

		// dot product gives cosine
		t_real dC = inner(vec0, vec1);
		t_real dAngle = std::atan2(dS, dC);

		// get signed angle
		if(pvec_norm)
		{
			// see if the cross product points along the direction
			// of the given normal
			if(inner(veccross, *pvec_norm) < t_real{0})
				dAngle = -dAngle;
		}

		return dAngle;
	}

	throw std::runtime_error("angle: only implemented for size == 2 and size == 3.");
}


/**
 * get the normal vector to a polygon
 */
template<class t_vec, template<class...> class t_cont = std::vector>
t_vec get_poly_normal(t_cont<t_vec>& vecPoly)
requires is_vec<t_vec>
{
	using namespace tl2_ops;
	using t_real = typename t_vec::value_type;

	// line from centre to vertex
	const t_vec vecCentre = mean(vecPoly);
	// face normal
	t_vec vecNormBest = zero<t_vec>(vecCentre.size());
	t_real tBestCross = t_real(0);

	// find non-collinear vectors
	for(std::size_t iVecPoly=1; iVecPoly<vecPoly.size(); ++iVecPoly)
	{
		t_vec vecNorm = cross<t_vec>({vecPoly[0]-vecCentre, vecPoly[1]-vecCentre});
		t_real tCross = norm(vecNorm);
		if(tCross > tBestCross)
		{
			tBestCross = tCross;
			vecNormBest = vecNorm;
		}
	}

	// nothing found
	if(vecNormBest.size() < vecCentre.size())
		return t_vec{};

	return vecNormBest;
}



/**
 * sort vertices in a convex polygon
 */
template<class t_vec, template<class...> class t_cont = std::vector>
void sort_poly_verts_norm(t_cont<t_vec>& vecPoly, const t_vec& _vecNorm)
requires is_vec<t_vec>
{
	using namespace tl2_ops;

	if(vecPoly.size() <= 1)
		return;

	// line from centre to vertex
	const t_vec vecCentre = mean(vecPoly);
	const t_vec vecNorm = _vecNorm / norm(_vecNorm);

	t_vec vec0 = vecPoly[0] - vecCentre;

	std::stable_sort(vecPoly.begin(), vecPoly.end(),
	[&vecCentre, &vec0, &vecNorm](const t_vec& vertex1, const t_vec& vertex2) -> bool
	{
		t_vec vec1 = vertex1 - vecCentre;
		t_vec vec2 = vertex2 - vecCentre;

		return angle(vec0, vec1, &vecNorm) < angle(vec0, vec2, &vecNorm);
	});
}


/**
 * sort vertices in a convex polygon using an absolute centre for determining the normal
 * @returns normal
 */
template<class t_vec, template<class...> class t_cont = std::vector>
t_vec sort_poly_verts(t_cont<t_vec>& vecPoly, const t_vec& vecAbsCentre,
	bool make_norm_perp_to_poly = false)
requires is_vec<t_vec>
{
	using namespace tl2_ops;

	if(vecPoly.size() <= 1)
		return t_vec{};

	// line from centre to vertex
	const t_vec vecCentre = mean(vecPoly);
	// face normal
	t_vec vecNorm = vecCentre - vecAbsCentre;

	sort_poly_verts_norm<t_vec, t_cont>(vecPoly, vecNorm);

	if(make_norm_perp_to_poly)
	{
		t_vec normal = get_poly_normal<t_vec, t_cont>(vecPoly);

		// do they point in the same direction?
		if(inner(normal, vecNorm) < 0.)
			normal = -normal;
		vecNorm = normal;
	}

	vecNorm /= norm(vecNorm);
	return vecNorm;
}


/**
 * sort vertices in a convex polygon determining normal
 * @returns normal
 */
template<class t_vec, template<class...> class t_cont = std::vector>
t_vec sort_poly_verts(t_cont<t_vec>& vecPoly)
requires is_vec<t_vec>
{
	using namespace tl2_ops;

	if(vecPoly.size() <= 1)
		return;

	t_vec vecNorm = get_poly_normal<t_vec, t_cont>(vecPoly);
	sort_poly_verts_norm<t_vec, t_cont>(vecPoly, vecNorm, true);

	vecNorm /= norm(vecNorm);
	return vecNorm;
}

}



#ifdef __TLIBS2_USE_QHULL__

namespace tl2_qh {

/**
 * calculates the convex hull
 * @see https://github.com/t-weber/misc/blob/master/geo/qhulltst.cpp
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_cont<t_vec>>, t_cont<t_vec>, t_cont<typename t_vec::value_type>>
get_convexhull(const t_cont<t_vec>& vecVerts)
requires tl2::is_vec<t_vec>
{
	using namespace tl2_ops;

	using t_real = typename t_vec::value_type;
	using t_real_qh = double;
	using t_facetlist_iter = typename orgQhull::QhullLinkedList<orgQhull::QhullFacet>::iterator;
	using t_vertexset_iter = typename orgQhull::QhullSet<orgQhull::QhullVertex>::iterator;


	// copy vertices
	int dim = vecVerts[0].size();
	std::size_t len = vecVerts.size()*dim;
	std::unique_ptr<t_real_qh[]> mem{new t_real_qh[len]};

	std::size_t i=0;
	for(const t_vec& vert : vecVerts)
	{
		for(int d=0; d<dim; ++d)
		{
			mem[i] = t_real_qh(vert[d]);
			++i;
		}
	}


	orgQhull::Qhull qhull{"tlibs2", dim, int(vecVerts.size()), mem.get(), "Qt"};
	orgQhull::QhullFacetList facets = qhull.facetList();

	t_cont<t_cont<t_vec>> vecPolys;
	t_cont<t_vec> vecNormals;
	t_cont<t_real> vecDists;
	vecPolys.reserve(facets.size());
	vecNormals.reserve(facets.size());
	vecDists.reserve(facets.size());
	//const t_vec vecCentre = tl2::mean(vecVerts);

	for(t_facetlist_iter iter=facets.begin(); iter!=facets.end(); ++iter)
	{
		// triangulable?
		if(iter->isUpperDelaunay())
			continue;

		orgQhull::QhullVertexSet vertices = iter->vertices();

		t_cont<t_vec> vecPoly;
		vecPoly.reserve(vertices.size());

		for(t_vertexset_iter iterVertex=vertices.begin(); iterVertex!=vertices.end(); ++iterVertex)
		{
			orgQhull::QhullPoint point = (*iterVertex).point();
			t_vec vecPoint(dim);
			for(int i=0; i<dim; ++i)
				vecPoint[i] = t_real(point[i]);

			vecPoly.emplace_back(std::move(vecPoint));
		}


		t_vec vecNormal; //= tl2::sort_poly_verts<t_vec, t_cont>(vecPoly, vecCentre, true);
		vecNormal.resize(dim);

		orgQhull::QhullHyperplane plane = iter->hyperplane();
		const t_real_qh* planenorm = plane.coordinates();
		const t_real_qh planedist = plane.offset();
		for(int i=0; i<dim; ++i)
			vecNormal[i] = t_real(planenorm[i]);

		vecPolys.emplace_back(std::move(vecPoly));
		vecNormals.emplace_back(std::move(vecNormal));
		vecDists.emplace_back(planedist);
	}

	// too few polygons => remove polyhedron
	if(vecPolys.size() < 3)
		vecPolys = decltype(vecPolys){};

	return std::make_tuple(vecPolys, vecNormals, vecDists);
}
}
#endif



namespace tl2 {
// ----------------------------------------------------------------------------
// Quaternions
// @see (Kuipers 2002) for infos
// ----------------------------------------------------------------------------

template<class t_quat> t_quat unit_quat() requires is_quat<t_quat>
{
	return t_quat(1, 0,0,0);
}


/**
 * calculates the quaternion inverse
 * @see (Bronstein 2008), Ch. 4
 * @see (Kuipers 2002), p. 112
 * @see https://en.wikipedia.org/wiki/Quaternion#Conjugation,_the_norm,_and_reciprocal
 */
template<class t_quat> t_quat inv(const t_quat& q) requires is_quat<t_quat>
{
	t_quat qc{q.R_component_1(), -q.R_component_2(), -q.R_component_3(), -q.R_component_4()};
	return qc / (q*qc);
}


/**
 * quaternion product
 * @see (Kuipers 2002), pp. 106-110
 * @see https://en.wikipedia.org/wiki/Quaternion#Scalar_and_vector_parts
 */
template<class t_quat, class t_vec>
t_quat prod(const t_quat& q1, const t_quat& q2)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using T = typename t_quat::value_type;

	T r1 = q1.R_component_1();
	T r2 = q2.R_component_1();

	t_vec vec1 = create<t_vec>({q1.R_component_2(), q1.R_component_3(), q1.R_component_4()});
	t_vec vec2 = create<t_vec>({q2.R_component_2(), q2.R_component_3(), q2.R_component_4()});

	T r = r1*r2 - inner<t_vec>(vec1, vec2);
	t_vec vec = r1*vec2 + r2*vec1 + cross<t_vec>({vec1, vec2});;

	return t_quat(r, vec[0], vec[1], vec[2]);
}


/**
 * 3x3 matrix -> quat
 * @desc algo from: http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q55
 */
template<class t_mat, class t_quat>
t_quat rot3_to_quat(const t_mat& rot)
requires is_quat<t_quat> && is_mat<t_mat>
{
	using T = typename t_quat::value_type;
	const T tr = trace<t_mat>(rot);
	T v[3], w;

	if(tr > T(0))	// scalar component is largest
	{
		w = T(0.5) * std::sqrt(tr+T(1));
		v[0] = (rot(2,1) - rot(1,2)) / (T(4)*w);
		v[1] = (rot(0,2) - rot(2,0)) / (T(4)*w);
		v[2] = (rot(1,0) - rot(0,1)) / (T(4)*w);
	}
	else
	{
		for(std::size_t iComp=0; iComp<3; ++iComp)	// find largest vector component
		{
			const std::size_t iM = iComp;		// major comp.
			const std::size_t im1 = (iComp+1)%3;	// minor comp. 1
			const std::size_t im2 = (iComp+2)%3;	// minor comp. 2

			if(rot(iM,iM)>=rot(im1,im1) && rot(iM,iM)>=rot(im2,im2))
			{
				v[iM] = T(0.5) * std::sqrt(T(1) + rot(iM,iM) - rot(im1,im1) - rot(im2,im2));
				v[im1] = (rot(im1, iM) + rot(iM, im1)) / (v[iM]*T(4));
				v[im2] = (rot(iM, im2) + rot(im2, iM)) / (v[iM]*T(4));
				w = (rot(im2,im1) - rot(im1,im2)) / (v[iM]*T(4));

				break;
			}

			if(iComp>=2)
				throw std::runtime_error("rot3_to_quat: invalid condition.");
		}
	}

	t_quat quatRet{w, v[0],v[1],v[2]};
	T norm_eucl = std::sqrt(w*w + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	return quatRet / norm_eucl;
}


/**
 * quat -> 3x3 matrix
 * @see (Desktop Bronstein 2008), formulas (4.162a/b)
 * @see (Bronstein 2008), p. 296, formulas (4.109a/b)
 */
template<class t_quat, class t_mat>
t_mat quat_to_rot3(const t_quat& quat)
requires is_quat<t_quat> && is_mat<t_mat>
{
	t_quat qc{quat.R_component_1(), -quat.R_component_2(), -quat.R_component_3(), -quat.R_component_4()};
	const t_quat i{0,1,0,0}, j{0,0,1,0}, k{0,0,0,1};
	const t_quat cols[] = { quat*i*qc, quat*j*qc, quat*k*qc };

	t_mat mat = unit<t_mat>(3);
	for(std::size_t icol=0; icol<3; ++icol)
	{
		mat(0, icol) = cols[icol].R_component_2();
		mat(1, icol) = cols[icol].R_component_3();
		mat(2, icol) = cols[icol].R_component_4();
	}

	return mat;
}


/**
 * vector -> quat
 * @see (Kuipers 2002), p. 114
 */
template<class t_vec, class t_quat>
t_quat vec3_to_quat(const t_vec& vec)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using T = typename t_vec::value_type;
	return t_quat{T(0), vec[0], vec[1], vec[2]};
}


/**
 * quat, vector product
 * @see (Kuipers 2002), p. 127
 */
template<class t_quat, class t_vec>
t_vec quat_vec_prod(const t_quat& q, const t_vec& v)
requires is_quat<t_quat> && is_vec<t_vec>
{
	t_quat qv = vec3_to_quat<t_vec, t_quat>(v);
	t_quat qc{q.R_component_1(), -q.R_component_2(), -q.R_component_3(), -q.R_component_4()};
	t_quat qvq =  q * qv * qc;

	return create<t_vec>({ qvq.R_component_2(), qvq.R_component_3(), qvq.R_component_4() });
}


/**
 * quat -> complex 2x2 matrix
 * @see (Scherer 2010), p.173
 * @see (DesktopBronstein08), ch. 4, equations (4.163a) and (4.163b)
 * @see (Bronstein08), ch. 4, p. 296, equation (4.110a) and (4.110b)
 */
template<class t_mat, class t_mat_cplx, class t_quat>
t_mat_cplx quat_to_cmat(const t_quat& quat)
requires is_quat<t_quat> && is_mat<t_mat> && is_mat<t_mat_cplx>
{
	using t_cplx = typename t_mat_cplx::value_type;

	const auto matI = unit<t_mat_cplx>(2);
	t_mat_cplx mat =
		t_cplx(quat.R_component_1()) * matI +
		t_cplx(quat.R_component_2()) * su2_matrix<t_mat_cplx>(0) +
		t_cplx(quat.R_component_3()) * su2_matrix<t_mat_cplx>(1) +
		t_cplx(quat.R_component_4()) * su2_matrix<t_mat_cplx>(2);

	return mat;
}


/**
 * rotation angle
 * @see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
 */
template<class t_quat>
typename t_quat::value_type rotation_angle(const t_quat& quat)
requires is_quat<t_quat>
{
	using t_real = typename t_quat::value_type;
	return t_real{2}*std::acos(quat.R_component_1());
}


/**
 * quat -> rotation axis
 * @see (Bronstein 2008), Ch. 4, pp. 301-302
 * @see https://en.wikipedia.org/wiki/Quaternion#Exponential,_logarithm,_and_power_functions
 * @see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
 */
template<class t_quat, class t_vec>
std::pair<t_vec, typename t_vec::value_type> rotation_axis(const t_quat& quat)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec vec = create<t_vec>({
		quat.R_component_2(),
		quat.R_component_3(),
		quat.R_component_4()
	});

	T angle = rotation_angle(quat);
	vec /= std::sin(T{0.5}*angle);

	return std::make_pair(vec, angle);
}


/**
 * rotation axis -> quat
 * @see (Desktop Bronstein 2008), formula (4.193)
 * @see (Bronstein 2008), ch. 4, pp. 301-302
 * @see https://en.wikipedia.org/wiki/Quaternion#Exponential,_logarithm,_and_power_functions
 * @see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
 */
template<class t_vec, class t_quat>
t_quat rotation_quat(const t_vec& vec, typename t_vec::value_type angle)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	const T s = std::sin(T(0.5)*angle);
	const T c = std::cos(T(0.5)*angle);
	const T n = norm<t_vec>(vec);

	const T x = s * vec[0] / n;
	const T y = s * vec[1] / n;
	const T z = s * vec[2] / n;
	const T r = c;

	return t_quat{r, x,y,z};
}


/**
 * quaternion to rotate vec0 into vec1
 */
template<class t_vec, class t_quat>
t_quat rotation_quat(const t_vec& _vec0, const t_vec& _vec1)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec vec0 = _vec0 / norm<t_vec>(_vec0);
	t_vec vec1 = _vec1 / norm<t_vec>(_vec1);

	// parallel vectors -> do nothing
	if(equals<t_vec>(vec0, vec1))
	{
		return unit_quat<t_quat>();
	}

	// antiparallel vectors -> rotate about any perpendicular axis
	else if(equals<t_vec>(vec0, -vec1))
	{
		t_vec vecPerp = perp<t_vec>(vec0);
		return rotation_quat<t_vec, t_quat>(vecPerp, pi<T>);
	}

	// rotation axis from cross product
	t_vec vecaxis = cross<t_vec>({vec0, vec1});

	T dC = inner<t_vec>(vec0, vec1);
	T dS = norm<t_vec>(vecaxis);

	// rotation angle
	T dAngle = std::atan2(dS, dC);

	return rotation_quat<t_vec, t_quat>(vecaxis, dAngle);
}


template<class t_quat>
t_quat rotation_quat_x(typename t_quat::value_type angle)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	return t_quat{std::cos(T(0.5)*angle), std::sin(T(0.5)*angle), T(0), T(0)};
}


template<class t_quat>
t_quat rotation_quat_y(typename t_quat::value_type angle)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	return t_quat{std::cos(T(0.5)*angle), T(0), std::sin(T(0.5)*angle), T(0)};
}


template<class t_quat>
t_quat rotation_quat_z(typename t_quat::value_type angle)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	return t_quat{std::cos(T(0.5)*angle), T(0), T(0), std::sin(T(0.5)*angle)};
}



/**
 * XYZ euler angles -> quat
 * @see (Kuipers 2002), pp. 166-167
 */
template<class t_quat>
t_quat euler_to_quat_xyz(
	typename t_quat::value_type phi, typename t_quat::value_type theta, typename t_quat::value_type psi)
requires is_quat<t_quat>
{
	t_quat q1 = rotation_quat_x<t_quat>(phi);
	t_quat q2 = rotation_quat_y<t_quat>(theta);
	t_quat q3 = rotation_quat_z<t_quat>(psi);

	return q3 * q2 * q1;
}


/**
 * ZXZ euler angles -> quat
 * @see (Kuipers 2002), pp. 166-167
 */
template<class t_quat>
t_quat euler_to_quat_zxz(
	typename t_quat::value_type phi, typename t_quat::value_type theta, typename t_quat::value_type psi)
requires is_quat<t_quat>
{
	t_quat q1 = rotation_quat_z<t_quat>(phi);
	t_quat q2 = rotation_quat_x<t_quat>(theta);
	t_quat q3 = rotation_quat_z<t_quat>(psi);

	return q3 * q2 * q1;
}


/**
 * quat -> XYZ euler angles
 * @see http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
 */
template<class t_quat>
std::vector<typename t_quat::value_type>
quat_to_euler_xyz(const t_quat& quat)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	T q[] = { quat.R_component_1(), quat.R_component_2(),
		quat.R_component_3(), quat.R_component_4() };

	// formulas from:
	// http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
	T phi = std::atan2(T(2)*(q[0]*q[1] + q[2]*q[3]), T(1)-T(2)*(q[1]*q[1] + q[2]*q[2]));
	T theta = std::asin(T(2)*(q[0]*q[2] - q[3]*q[1]));
	T psi = std::atan2(T(2)*(q[0]*q[3] + q[1]*q[2]), T(1)-T(2)*(q[2]*q[2] + q[3]*q[3]));

	return std::vector<T>({ phi, theta, psi });
}


/**
 * @see (Desktop Bronstein 2008), formula (4.217)
 * @see (Bronstein 2008), p. 308, formula (4.165)
 */
template<class t_quat>
t_quat stereo_proj(const t_quat& quat)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	return (T{1}+quat) / (T{1}-quat);
}


/**
 * @see (Desktop Bronstein 2008), formula (4.217)
 * @see (Bronstein 2008), p. 308, formula (4.165)
 */
template<class t_quat>
t_quat stereo_proj_inv(const t_quat& quat)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	return (T{1}-quat) / (T{1}+quat);
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Fourier transform
// @see (Scarpino 2011), ch. 14 for infos.
// ----------------------------------------------------------------------------

/**
 * dft
 * @see http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
 */
template<typename T = double, class t_cplx = std::complex<T>,
	template<class...> class t_cont = std::vector>
requires is_complex<t_cplx>
t_cplx dft_coeff(int k, const t_cont<t_cplx>& invec, bool bInv = false)
{
	const std::size_t N = invec.size();

	t_cplx imag(0., 1.);
	t_cplx f(0., 0.);

	for(std::size_t j=0; j<N; ++j)
	{
		T dv = T(-2)*pi<T>*T(j)*T(k)/T(N);
		if(bInv) dv = -dv;
		f += invec[j] * (std::cos(dv) + imag*std::sin(dv));
	}

	return f;
}


/**
 * dft
 * @see http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
 */
template<typename T = double, class t_cplx = std::complex<T>,
	template<class...> class t_cont = std::vector>
requires is_complex<t_cplx>
t_cont<t_cplx> dft(const t_cont<t_cplx>& invec,
	bool bInv = false, bool bNorm = false)
{
	const std::size_t N = invec.size();
	t_cont<t_cplx> outvec;
	outvec.resize(N);

	for(std::size_t k=0; k<N; ++k)
	{
		outvec[k] = dft_coeff<T, t_cplx, t_cont>(k, invec, bInv);

		if(bNorm && bInv)
			outvec[k] /= N;
	}

	return outvec;
}


/**
 * fft
 * @see (Scarpino 2011), ch. 14.
 */
template<typename T = double, class t_cplx = std::complex<T>>
requires is_complex<t_cplx>
t_cplx fft_factor(T N, T k, bool bInv = false)
{
	T ph = bInv ? -1 : 1.;

	T c = std::cos(T(2)*pi<T>*k/N * ph);
	T s = std::sin(T(2)*pi<T>*k/N * ph);

	return t_cplx{c, -s};
}


/**
 * fft
 * @see (Scarpino 2011), ch. 14.
 */
template<typename T = double, class t_cplx = std::complex<T>,
	template<class...> class t_cont = std::vector>
requires is_complex<t_cplx>
t_cont<t_cplx> fft_reorder(const t_cont<t_cplx>& vecIn)
{
	t_cont<std::size_t> vecIdx =
		bit_reverse_indices<std::size_t, t_cont>(vecIn.size());

	t_cont<t_cplx> vecInRev;
	vecInRev.reserve(vecIn.size());

	for(std::size_t i=0; i<vecIn.size(); ++i)
		vecInRev.push_back(vecIn[vecIdx[i]]);

	return vecInRev;
}


/**
 * fft
 * @see (Scarpino 2011), ch. 14.
 */
template<typename T = double, class t_cplx = std::complex<T>,
	template<class...> class t_cont = std::vector>
requires is_complex<t_cplx>
t_cont<t_cplx> fft_merge(const t_cont<t_cplx>& vecIn, bool bInv = false)
{
	const std::size_t N = vecIn.size();
	const std::size_t N2 = N/2;

	if(N==0 || N==1)
		return vecIn;

	auto split_vec = [](const t_cont<t_cplx>& vec)
		-> std::pair<t_cont<t_cplx>, t_cont<t_cplx>>
	{
		std::size_t N = vec.size();

		t_cont<t_cplx> vec1, vec2;
		vec1.reserve(N/2);
		vec2.reserve(N/2);

		for(std::size_t i=0; i<N/2; ++i)
		{
			vec1.push_back(vec[i]);
			vec2.push_back(vec[N/2 + i]);
		}

		return std::make_pair(std::move(vec1), std::move(vec2));
	};

	auto pair = split_vec(vecIn);
	t_cont<t_cplx> vec1 = fft_merge<T, t_cplx, t_cont>(pair.first, bInv);
	t_cont<t_cplx> vec2 = fft_merge<T, t_cplx, t_cont>(pair.second, bInv);

	t_cont<t_cplx> vecOut;
	vecOut.resize(N);

	for(std::size_t i=0; i<N2; ++i)
	{
		vecOut[i] = vec1[i] + vec2[i]*fft_factor<T, t_cplx>(N, i, bInv);
		vecOut[N2+i] = vec1[i] + vec2[i]*fft_factor<T, t_cplx>(N, N2+i, bInv);
	}

	return vecOut;
}


/**
 * fft
 * @see (Scarpino 2011), ch. 14.
 */
template<typename T = double, class t_cplx = std::complex<T>,
	template<class...> class t_cont = std::vector>
requires is_complex<t_cplx>
t_cont<t_cplx> fft(const t_cont<t_cplx>& vecIn,
	bool bInv = false, bool bNorm = false)
{
	const std::size_t n = vecIn.size();
	t_cont<t_cplx> vecOut;
	vecOut.resize(n);

	vecOut = fft_reorder<T, t_cplx, t_cont>(vecIn);
	vecOut = fft_merge<T, t_cplx, t_cont>(vecOut, bInv);

	if(bInv && bNorm)
		for(t_cplx& c : vecOut)
			c /= n;

	return vecOut;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// numerical integration and differentiation
// ----------------------------------------------------------------------------
/**
 * numerical differentiation
 */
template<class t_cont_out = std::vector<double>, class t_cont_in = std::vector<double>>
t_cont_out diff(const t_cont_in& xs, const t_cont_in& ys)
{
	using t_real = typename t_cont_in::value_type;

	const std::size_t N = xs.size();
	t_cont_out new_ys{};
	new_ys.reserve(N);

	for(std::size_t i=0; i<N-1; ++i)
	{
		t_real diff_val = (ys[i+1]-ys[i]) / (xs[i+1]-xs[i]);
		new_ys.push_back(diff_val);
	}

	// repeat last value
	new_ys[N-1] = new_ys[N-2];
	return new_ys;
}


/**
 * trapezoid rule
 * @see https://en.wikipedia.org/wiki/Trapezoidal_rule
 */
template<class R=double, class A=double>
R numint_trap(const std::function<R(A)>& fkt,
	A x0, A x1)
{
	return R(0.5)*R(x1-x0) * (fkt(x0) + fkt(x1));
}


/**
 * trapezoid rule
 * @see https://en.wikipedia.org/wiki/Trapezoidal_rule
 */
template<class R=double, class A=double>
R numint_trapN(const std::function<R(A)>& fkt,
	A x0, A x1, std::size_t N)
{
	const A xstep = A(x1-x0)/A(N);

	R xsum = fkt(x0) + fkt(x1);
	for(std::size_t i=1; i<N; ++i)
		xsum += R(2)*fkt(x0 + A(i)*xstep);

	xsum *= R(0.5)*R(xstep);
	return R(xsum);
}


/**
 * rectangle rule
 * @see https://en.wikipedia.org/wiki/Rectangle_method
 */
template<class R=double, class A=double>
R numint_rect(const std::function<R(A)>& fkt,
	A x0, A x1, std::size_t N)
{
	const A xstep = (x1-x0)/A(N);

	R xsum = R(0);
	for(std::size_t i=0; i<N; ++i)
		xsum += fkt(x0 + A(i)*xstep);

	xsum *= R(xstep);
	return xsum;
}


/**
 * Simpson's rule
 * @see https://en.wikipedia.org/wiki/Simpson%27s_rule
 */
template<class R=double, class A=double>
R numint_simp(const std::function<R(A)>& fkt, A x0, A x1)
{
	return (fkt(x0) + 4.*fkt(0.5*(x0+x1)) + fkt(x1)) * (x1-x0)/6.;
}


/**
 * Simpson's rule
 * @see https://en.wikipedia.org/wiki/Simpson%27s_rule
 */
template<class R=double, class A=double>
R numint_simpN(const std::function<R(A)>& fkt, A x0, A x1, std::size_t N)
{
	const A xstep = (x1-x0)/A(N);
	R xsum = fkt(x0) + fkt(x1);

	for(std::size_t i=1; i<=N/2; ++i)
	{
		xsum += R(2) * fkt(x0 + A(2*i)*xstep);
		xsum += R(4) * fkt(x0 + A(2*i-1)*xstep);
	}
	xsum -= R(2)*fkt(x0 + A(2*N/2)*xstep);

	xsum *= R(xstep)/R(3);
	return xsum;
}


/**
 * convolution integral of fkt0 and fkt1
 * @see https://en.wikipedia.org/wiki/Convolution
 */
template<class R=double, class A=double>
R convolute(const std::function<R(A)>& fkt0,
	const std::function<R(A)>& fkt1,
	A x, A x0, A x1, std::size_t N)
{
	std::function<R(A,A)> fkt = [&fkt0, &fkt1](A t, A tau) -> R
	{
		return fkt0(tau) * fkt1(t-tau);
	};

	// ... at fixed arg x
	std::function<R(A)> fktbnd = [&fkt, &x](A t) -> R { return fkt(x, t); };

	return numint_simpN(fktbnd, x0, x1, N);
}


/**
 * discrete convolution
 * @see https://en.wikipedia.org/wiki/Convolution#Discrete_convolution
 */
template<class cont_type = std::vector<double>>
cont_type convolute_discrete(const cont_type& f, const cont_type& g)
{
	const std::size_t M = f.size();
	const std::size_t N = g.size();

	cont_type conv;
	conv.reserve(M+N-1);

	for(std::size_t n=0; n<M+N-1; ++n)
	{
		typename cont_type::value_type val = 0.;

		for(std::size_t m=0; m<M; ++m)
			if(n>=m && n-m<N)
				val += f[m] * g[n-m];

		conv.push_back(val);
	}

	return conv;
}


/**
 * Newton iteration
 * @see https://en.wikipedia.org/wiki/Newton%27s_method
 */
template<class T = double>
T newton(const std::function<T(T)>& fkt, const std::function<T(T)>& diff,
	T x, std::size_t imax = 128, T eps = std::numeric_limits<T>::epsilon())
{
	T xnew = x;
	for(std::size_t i=0; i<imax; ++i)
	{
		xnew = x - fkt(x)/diff(x);

		if(std::abs(xnew - x) < eps)
			break;

		x = xnew;
	}

	return xnew;
}
// ----------------------------------------------------------------------------


}

// ----------------------------------------------------------------------------
#endif
