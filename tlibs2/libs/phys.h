/**
 * tlibs2
 * physics library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2012-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * @note Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @note Additional functions were forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 * @note Further functions and updates forked on 1-Feb-2021 and 19-Apr-2021 from my privately developed "geo" and "misc" projects (https://github.com/t-weber/geo and https://github.com/t-weber/misc).
 *
 * @note for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "magtools", "geo", and "misc" projects
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#ifndef __TLIBS2_PHYS__
#define __TLIBS2_PHYS__

#include "units.h"
#include "log.h"
#include "maths.h"

#include <boost/units/pow.hpp>


namespace tl2 {


// --------------------------------------------------------------------------------
// constants
// --------------------------------------------------------------------------------
// import scipy.constants as co
// E2KSQ = 2.*co.neutron_mass/(co.Planck/co.elementary_charge*1000./2./co.pi)**2. / co.elementary_charge*1000. * 1e-20
template<class T=double> constexpr T KSQ2E = T(0.5) * hbar<T>/angstrom<T>/m_n<T> * hbar<T>/angstrom<T>/meV<T>;
template<class T=double> constexpr T E2KSQ = T(1) / KSQ2E<T>;
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------

/**
 * differentiated Bragg equation:
 * n lam = 2d sin(th)				| diff
 * n dlam = 2dd sin(th) + 2d cos(th) dth	| / Bragg equ
 * dlam/lam = dd/d + cos(th)/sin(th) dth
 *
 * n G = 2k sin(th)
 * n dG = 2dk sin(th) + 2k cos(th) dth
 * dG/G = dk/k + cos(th)/sin(th) dth
 */
template<class Sys, class Y=double>
Y bragg_diff(Y dDoverD, const t_angle<Sys,Y>& theta, Y dTheta)
{
	Y dLamOverLam = dDoverD + units::cos(theta)/units::sin(theta) * dTheta;
	return dLamOverLam;
}


/**
 * kinematic plane
 * @see (ILL Neutron Data Booklet), sec. 2.6-2
 *
 * Q_vec = ki_vec - kf_vec
 * Q^2 = ki^2 + kf^2 - 2ki kf cos 2th	| * hbar^2 / (2 mn)
 *
 * using: Ei = hbar^2 ki^2 / (2 mn)
 * Q^2 * hbar^2 / (2 mn) = Ei + Ef - 2 ki kf cos(2th) * hbar^2 / (2 mn)
 *
 * using: ki^2 = 2 mn Ei / hbar^2
 * Q^2 = [Ei + Ef - 2 sqrt(Ei) sqrt(Ef) cos 2th] * 2 mn / hbar^2
 *
 * using: dE = Ei - Ef, Ef = Ei - dE
 * Q^2 = [2 Ei - dE - 2 sqrt(Ei (Ei - dE)) cos 2th] * 2 mn / hbar^2
 */
template<class Sys, class Y>
t_wavenumber<Sys,Y> kinematic_plane(bool bFixedKi,
	const t_energy<Sys,Y>& EiEf, const t_energy<Sys,Y>& DeltaE,
	const t_angle<Sys,Y>& twotheta)
{
	t_energy<Sys,Y> dE = DeltaE;
	if(bFixedKi)
		dE = -dE;

	auto c = Y(2.)*m_n<Y> / (hbar<Y>*hbar<Y>);
	t_wavenumber<Sys,Y> Q =
		units::sqrt(c *
		(Y(2.)*EiEf + dE - Y(2.)*units::cos(twotheta) *
		units::sqrt(EiEf*(EiEf + dE))));

	return Q;
}


/**
 * kinematic plane
 * @see (ILL Neutron Data Booklet), sec. 2.6-2
 *
 * solving the above equation for dE using sage:
 *   Q, Ei, dE, ctt, c = var("Q, Ei, dE, ctt, c")
 *   equ = (Q^2 -2*Ei*c + dE*c)^2 == 4*Ei*(Ei-dE)*c^2*ctt^2
 *   equ.solve(dE)
 */
template<class Sys, class Y>
t_energy<Sys,Y> kinematic_plane(bool bFixedKi, bool bBranch,
	const t_energy<Sys,Y>& EiEf, const t_wavenumber<Sys,Y>& Q,
	const t_angle<Sys,Y>& twotheta)
{
	auto c = Y(2.)*m_n<Y> / (hbar<Y>*hbar<Y>);
	auto c2 = c*c;

	auto EiEf2 = EiEf*EiEf;

	Y ctt = units::cos(twotheta);
	Y ctt2 = ctt*ctt;

	Y dSign = bBranch ? Y(1.) : Y(-1.);
	Y dSignFixedKf = bFixedKi ? Y(-1.) : Y(1.);

	t_energy<Sys,Y> dE =
			dSignFixedKf*Y(2.) * EiEf * ctt2
			- dSignFixedKf*Y(2.) * EiEf
			+ dSignFixedKf * Q*Q / c
			+ dSign*Y(2.) * ctt/c * units::sqrt(c2*ctt2*EiEf2 - c2*EiEf2 + c*EiEf*Q*Q);

	return dE;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// tas calculations
// @see M. D. Lumsden, J. L. Robertson, and M. Yethiraj, J. Appl. Crystallogr. 38(3), pp. 405–411 (2005), doi: 10.1107/S0021889805004875.
// @see (Shirane 2002), Ch. 1.3
// --------------------------------------------------------------------------------

/**
 * angle between ki and kf in the scattering triangle
 * @returns nullopt if the angle can't be reached
 *
 * |Q> = |ki> - |kf>
 * Q^2 = ki^2 + kf^2 - 2*<ki|kf>
 * 2*<ki|kf> = ki^2 + kf^2 - Q^2
 * cos phi = (ki^2 + kf^2 - Q^2) / (2 ki*kf)
 */
template<typename t_real>
std::optional<t_real> calc_tas_angle_ki_kf(
	t_real ki, t_real kf, t_real Q, t_real sense=1)
{
	t_real c = (ki*ki + kf*kf - Q*Q) / (t_real(2)*ki*kf);
	if(std::abs(c) > t_real(1))
		return std::nullopt;
	return sense*std::acos(c);
}


/**
 * angle between ki and Q in the scattering triangle
 * @returns nullopt if the angle can't be reached
 *
 * |Q> = |ki> - |kf>
 * |kf> = |ki> + |Q>
 * kf^2 = ki^2 + Q^2 - 2*<ki|Q>
 * 2*<ki|Q> = ki^2 + Q^2 - kf^2
 * cos phi = (ki^2 + Q^2 - kf^2) / (2 ki*Q)
 */
template<typename t_real>
std::optional<t_real> calc_tas_angle_ki_Q(
	t_real ki, t_real kf, t_real Q, t_real sense=1)
{
	t_real c = (ki*ki + Q*Q - kf*kf) / (t_real(2)*ki*Q);
	if(std::abs(c) > t_real(1))
		return std::nullopt;
	return sense*std::acos(c);
}


/**
 * angle between ki and Q in the scattering triangle
 * (version with units)
 */
template<class Sys, class Y>
t_angle<Sys,Y> calc_tas_angle_ki_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf,
	const t_wavenumber<Sys,Y>& Q,
	bool bPosSense=1, bool bAngleOutsideTriag=0)
{
	t_angle<Sys,Y> angle;

	if(Q*angstrom<Y> == Y(0.))
	{
		angle = pi<Y>/Y(2) * radians<Y>;
	}
	else
	{
		auto c = (ki*ki - kf*kf + Q*Q) / (Y(2.)*ki*Q);
		if(units::abs(c) > Y(1.))
			throw Err("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(bAngleOutsideTriag) angle = pi<Y>*radians<Y> - angle;
	if(!bPosSense) angle = -angle;

	return angle;
}


/**
 * angle between kf and Q in the scattering triangle
 * (version with units)
 *
 * Q_vec = ki_vec - kf_vec
 * ki_vec = Q_vec + kf_vec
 * ki^2 = Q^2 + kf^2 + 2Q kf cos th
 * cos th = (ki^2 - Q^2 - kf^2) / (2Q kf)
 */
template<class Sys, class Y>
t_angle<Sys,Y> calc_tas_angle_kf_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf,
	const t_wavenumber<Sys,Y>& Q,
	bool bPosSense=1, bool bAngleOutsideTriag=1)
{
	t_angle<Sys,Y> angle;

	if(Q*angstrom<Y> == Y(0.))
		angle = pi<Y>/Y(2) * radians<Y>;
	else
	{
		auto c = (ki*ki - kf*kf - Q*Q) / (Y(2.)*kf*Q);
		if(units::abs(c) > Y(1.))
			throw Err("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(!bAngleOutsideTriag) angle = pi<Y>*radians<Y> - angle;
	if(!bPosSense) angle = -angle;

	return angle;
}


/**
 * get length of Q
 * |Q> = |ki> - |kf>
 * Q^2 = ki^2 + kf^2 - 2*<ki|kf>
 * Q^2 = ki^2 + kf^2 - 2*ki*kf*cos(a4)
 */
template<typename t_real>
t_real calc_tas_Q_len(t_real ki, t_real kf, t_real a4)
{
	t_real Qsq = ki*ki + kf*kf - t_real(2)*ki*kf*std::cos(a4);
	return std::sqrt(Qsq);
}


/**
 * get length of Q
 * (version with units)
 * |Q> = |ki> - |kf>
 * Q^2 = ki^2 + kf^2 - 2*<ki|kf>
 * Q^2 = ki^2 + kf^2 - 2*ki*kf*cos(a4)
 */
template<class Sys, class Y>
t_wavenumber<Sys,Y>
calc_tas_Q_len(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf, const t_angle<Sys,Y>& tt)
{
	t_dimensionless<Sys,Y> ctt = units::cos(tt);
	decltype(ki*ki) Qsq = ki*ki + kf*kf - Y(2.)*ki*kf*ctt;

	if(Y(Qsq*angstrom<Y>*angstrom<Y>) < Y(0.))
	{
		// TODO
		Qsq = -Qsq;
	}

	//t_wavenumber<Sys,Y> Q = units::sqrt(Qsq);
	t_wavenumber<Sys,Y> Q = my_units_sqrt<t_wavenumber<Sys,Y>>(Qsq);
	return Q;
}


/**
 * get tas a4 angles
 * (version with units)
 * Q_vec = ki_vec - kf_vec
 * Q^2 = ki^2 + kf^2 - 2ki kf cos 2th
 * cos 2th = (-Q^2 + ki^2 + kf^2) / (2ki kf)
 */
template<class Sys, class Y>
t_angle<Sys,Y> calc_tas_a4(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf, const t_wavenumber<Sys,Y>& Q,
	bool bPosSense=1)
{
	t_dimensionless<Sys,Y> ttCos = (ki*ki + kf*kf - Q*Q)/(Y(2.)*ki*kf);
	if(units::abs(ttCos) > Y(1.))
		throw Err("Scattering triangle not closed.");

	t_angle<Sys,Y> tt;
	tt = units::acos(ttCos);

	if(!bPosSense) tt = -tt;
	return tt;
}


/**
 * get tas a3 and a4 angles
 * @return [a3, a4, distance of Q to the scattering plane]
 * @see M. D. Lumsden, et al., doi: 10.1107/S0021889805004875.
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
std::tuple<bool, t_real, t_real, t_real> calc_tas_a3a4(
	const t_mat& B, t_real ki_lab, t_real kf_lab,
	const t_vec& Q_rlu, const t_vec& orient_rlu, const t_vec& orient_up_rlu,
	t_real sample_sense = 1, t_real a3_offs = pi<t_real>)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	// metric from crystal B matrix
	t_mat G = tl2::metric<t_mat>(B);

	// length of Q vector
	t_real Q_len_lab = norm<t_mat, t_vec>(G, Q_rlu);

	// angle xi between Q and orientation reflex
	t_real xi = angle<t_mat, t_vec>(G, Q_rlu, orient_rlu);

	// sign/direction of xi
	t_vec xivec = cross<t_mat, t_vec>(G, orient_rlu, Q_rlu);
	t_real xidir = inner<t_mat, t_vec>(G, xivec, orient_up_rlu);
	if(xidir < t_real(0))
		xi = -xi;

	// angle psi between ki and Q
	std::optional<t_real> psi =
		calc_tas_angle_ki_Q<t_real>(ki_lab, kf_lab, Q_len_lab, sample_sense);
	if(!psi)
		return std::make_tuple(false, 0, 0, 0);

	// crystal and scattering angle
	t_real a3 = - *psi - xi + a3_offs;
	std::optional<t_real> a4 =
		calc_tas_angle_ki_kf<t_real>(ki_lab, kf_lab, Q_len_lab);
	if(!a4)
		return std::make_tuple(false, a3, 0, 0);
	*a4 *= sample_sense;

	// distance of Q to the scattering plane
	t_real dist_Q_plane = inner<t_mat, t_vec>(G, Q_rlu, orient_up_rlu);
	dist_Q_plane /= norm<t_mat, t_vec>(G, orient_up_rlu);

	return std::make_tuple(true, a3, *a4, dist_Q_plane);
}


/**
 * get hkl position of a tas
 * @return Q_rlu
 * @see M. D. Lumsden, et al., doi: 10.1107/S0021889805004875.
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
std::optional<t_vec> calc_tas_hkl(
	const t_mat& B, t_real ki_lab, t_real kf_lab, t_real Q_len_lab, t_real a3,
	const t_vec& orient_rlu, const t_vec& orient_up_rlu,
	t_real sample_sense = 1, t_real a3_offs = pi<t_real>)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	auto [Binv, ok] = inv<t_mat>(B);
	if(!ok)
		return std::nullopt;

	// angle psi between ki and Q
	std::optional<t_real> psi =
		calc_tas_angle_ki_Q<t_real>(ki_lab, kf_lab, Q_len_lab, sample_sense);
	if(!psi)
		return std::nullopt;

	// angle xi between Q and orientation reflex
	t_real xi = a3_offs - a3 - *psi;

	t_vec rotaxis_lab = B * orient_up_rlu;
	t_mat rotmat = rotation<t_mat, t_vec>(rotaxis_lab, xi, false);

	t_vec orient_lab = B * orient_rlu;
	t_vec Q_lab = rotmat * orient_lab;
	Q_lab /= norm<t_vec>(Q_lab);
	Q_lab *= Q_len_lab;

	t_vec Q_rlu = Binv * Q_lab;
	return Q_rlu;
}


/**
 * get a1 or a5 angle
 * @returns nullopt of the angle can't be reached
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 *
 * Bragg: n lam = 2d sin(theta)
 * n 2pi / k = 2d sin(theta)
 * n pi / k = d sin(theta)
 * theta = asin(n pi / (k d))
 */
template<class t_real>
std::optional<t_real> calc_tas_a1(t_real k, t_real d)
{
	t_real sintheta = pi<t_real> / (k*d);
	if(std::abs(sintheta) > t_real(1))
		return std::nullopt;
	return std::asin(sintheta);
}


/**
 * get a2 or a6 angle
 * (version with units)
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 */
template<class Sys, class Y>
t_angle<Sys,Y> calc_tas_a1(const t_wavenumber<Sys,Y>& k,
	const t_length<Sys,Y>& d, bool bPosSense=1)
{
	const Y order = Y(1.);
	t_length<Sys,Y> lam = Y(2.)*pi<Y> / k;
	auto dS = order*lam/(Y(2.)*d);
	if(std::abs(Y(dS)) > Y(1))
		throw Err("Invalid twotheta angle.");

	t_angle<Sys,Y> theta = units::asin(dS);
	if(!bPosSense)
		theta = -theta;
	return theta;
}


/**
 * get k from crystal angle
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 *
 * k = n pi / (d sin(theta))
 */
template<class t_real>
t_real calc_tas_k(t_real theta, t_real d)
{
	t_real sintheta = std::abs(std::sin(theta));
	return pi<t_real> / (d * sintheta);
}


/**
 * get k from crystal angle
 * (version with units)
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 */
template<class Sys, class Y>
t_wavenumber<Sys,Y> calc_tas_k(const t_angle<Sys,Y>& _theta,
	const t_length<Sys,Y>& d, bool bPosSense=1)
{
	t_angle<Sys,Y> theta = _theta;
	if(!bPosSense)
		theta = -theta;

	const Y order = Y(1.);

	// https://en.wikipedia.org/wiki/Bragg%27s_law
	t_length<Sys, Y> lam = Y(2.)*d/order * units::sin(theta);
	t_wavenumber<Sys,Y> k = Y(2.)*pi<Y> / lam;

	return k;
}


/**
 * get ki from kf and energy transfer
 */
template<class t_real>
t_real calc_tas_ki(t_real kf, t_real E)
{
	return std::sqrt(kf*kf + E2KSQ<t_real>*E);
}


/**
 * get kf from ki and energy transfer
 */
template<class t_real>
t_real calc_tas_kf(t_real ki, t_real E)
{
	return std::sqrt(ki*ki - E2KSQ<t_real>*E);
}


/**
 * get energy transfer from ki and kf
 */
template<class t_real>
t_real calc_tas_E(t_real ki, t_real kf)
{
	return (ki*ki - kf*kf) / E2KSQ<t_real>;
}


template<class Sys, class Y>
t_energy<Sys,Y> k2E(const t_wavenumber<Sys,Y>& k)
{
	Y dk = k*angstrom<Y>;
	Y dE = KSQ2E<Y> * dk*dk;
	return dE * meV<Y>;
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> E2k(const t_energy<Sys,Y>& _E, bool &bImag)
{
	bImag = (_E < Y(0.)*meV<Y>);
	t_energy<Sys,Y> E = bImag ? -_E : _E;
	const Y dE = E / meV<Y>;
	const Y dk = std::sqrt(E2KSQ<Y> * dE);
	return dk / angstrom<Y>;
}


/**
 * get energy transfer from ki and kf
 * (version with units)
 */
template<class Sys, class Y>
t_energy<Sys,Y> get_energy_transfer(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf)
{
	return k2E<Sys,Y>(ki) - k2E<Sys,Y>(kf);
}


/**
 * (hbar*ki)^2 / (2*mn)  -  (hbar*kf)^2 / (2mn)  =  E
 * 1) ki^2  =  +E * 2*mn / hbar^2  +  kf^2
 * 2) kf^2  =  -E * 2*mn / hbar^2  +  ki^2
 */
template<class Sys, class Y>
t_wavenumber<Sys,Y> get_other_k(const t_energy<Sys,Y>& E,
	const t_wavenumber<Sys,Y>& kfix, bool bFixedKi)
{
	auto kE_sq = E*Y(2.)*(m_n<Y>/hbar<Y>)/hbar<Y>;
	if(bFixedKi) kE_sq = -kE_sq;

	auto k_sq = kE_sq + kfix*kfix;
	if(k_sq*angstrom<Y>*angstrom<Y> < Y(0.))
		throw Err("Scattering triangle not closed.");

	//return units::sqrt(k_sq);
	return my_units_sqrt<t_wavenumber<Sys,Y>>(k_sq);
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------

/**
 * kf^3 mono/ana reflectivity factor
 * @see (Shirane 2002) p. 125
 */
template<class Sys, class Y>
Y ana_effic_factor(const t_wavenumber<Sys, Y>& kf, const t_angle<Sys, Y>& theta)
{
	return kf*kf*kf / units::tan(theta) * angstrom<Y>*angstrom<Y>*angstrom<Y>;
}


/**
 * kf^3 mono/ana reflectivity factor,
 * @see (Shirane 2002) p. 125
 */
template<class Sys, class Y>
Y ana_effic_factor(const t_wavenumber<Sys, Y>& kf, const t_length<Sys, Y>& d)
{
	t_angle<Sys, Y> theta = units::abs(calc_tas_a1<Sys, Y>(kf, d, true));
	return ana_effic_factor<Sys, Y>(kf, theta);
}

// --------------------------------------------------------------------------------



/**
 * Bose distribution
 * @see (Shirane 2002), p. 28
 * @see https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics
 */
template<class t_real=double>
t_real bose(t_real E, t_real T)
{
	const t_real _kB = kB<t_real> * kelvin<t_real>/meV<t_real>;

	t_real n = t_real(1)/(std::exp(std::abs(E)/(_kB*T)) - t_real(1));
	if(E >= t_real(0))
		n += t_real(1);

	return n;
}


/**
 * Bose factor with a lower cutoff energy
 * @see https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics
 */
template<class t_real=double>
t_real bose_cutoff(t_real E, t_real T, t_real E_cutoff=t_real(0.02))
{
	t_real dB;

	E_cutoff = std::abs(E_cutoff);
	if(std::abs(E) < E_cutoff)
		dB = bose<t_real>(sign(E)*E_cutoff, T);
	else
		dB = bose<t_real>(E, T);

	return dB;
}


/**
 * Bose factor
 * @see https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics
 */
template<class Sys, class Y>
Y bose(const t_energy<Sys,Y>& E, const t_temperature<Sys,Y>& T,
	t_energy<Sys,Y> E_cutoff = -meV<Y>)
{
	if(E_cutoff < Y(0)*meV<Y>)
		return bose<Y>(Y(E/meV<Y>), Y(T/kelvin<Y>));
	else
		return bose_cutoff<Y>(Y(E/meV<Y>), Y(T/kelvin<Y>),
			Y(E_cutoff/meV<Y>));
}


/**
 * DHO
 * @see B. Fak, B. Dorner, Physica B 234-236 (1997) pp. 1107-1108, doi: https://doi.org/10.1016/S0921-4526(97)00121-X
 */
template<class t_real=double>
t_real DHO_model(t_real E, t_real T, t_real E0, t_real hwhm, t_real amp, t_real offs)
{
	//if(E0*E0 - hwhm*hwhm < 0.) return 0.;
	return std::abs(bose<t_real>(E, T)*amp/(E0*pi<t_real>) *
		(hwhm/((E-E0)*(E-E0) + hwhm*hwhm) - hwhm/((E+E0)*(E+E0) + hwhm*hwhm)))
		+ offs;
}


// --------------------------------------------------------------------------------

/**
 * Fermi distribution
 * @see https://en.wikipedia.org/wiki/Fermi%E2%80%93Dirac_statistics
 */
template<class t_real=double>
t_real fermi(t_real E, t_real mu, t_real T)
{
	const t_real _kB = kB<t_real> * kelvin<t_real>/meV<t_real>;
	t_real n = t_real(1)/(std::exp((E-mu)/(_kB*T)) + t_real(1));
	return n;
}


/**
 * Fermi distribution
 * @see https://en.wikipedia.org/wiki/Fermi%E2%80%93Dirac_statistics
 */
template<class Sys, class Y>
Y fermi(const t_energy<Sys,Y>& E, const t_energy<Sys,Y>& mu,
	const t_temperature<Sys,Y>& T)
{
	return fermi<Y>(Y(E/meV<Y>), Y(mu/meV<Y>), Y(T/kelvin<Y>));
}

// --------------------------------------------------------------------------------


/**
 * get macroscopic from microscopic cross-section
 */
template<class Sys, class Y=double>
t_length_inverse<Sys, Y> macro_xsect(const t_area<Sys, Y>& xsect,
	unsigned int iNumAtoms, const t_volume<Sys, Y>& volUC)
{
	return xsect * Y(iNumAtoms) / volUC;
}



// --------------------------------------------------------------------------------

/**
 * thin lens equation: 1/f = 1/lenB + 1/lenA
 * @see https://en.wikipedia.org/wiki/Thin_lens
 */
template<class Sys, class Y=double>
t_length<Sys, Y> focal_len(const t_length<Sys, Y>& lenBefore, const t_length<Sys, Y>& lenAfter)
{
	const t_length_inverse<Sys, Y> f_inv = Y(1)/lenBefore + Y(1)/lenAfter;
	return Y(1) / f_inv;
}


/**
 * optimal mono/ana curvature,
 * @see (Shirane 2002) p. 66
 * @see NICOS: https://forge.frm2.tum.de/cgit/cgit.cgi/frm2/nicos/nicos-core.git/plain/nicos/devices/tas/mono.py
 * @see McStas: https://github.com/McStasMcXtrace/McCode/blob/master/mcstas-comps/optics/Monochromator_curved.comp
 */
template<class Sys, class Y=double>
t_length<Sys, Y> foc_curv(const t_length<Sys, Y>& lenBefore, const t_length<Sys, Y>& lenAfter,
	const t_angle<Sys, Y>& tt, bool bVert)
{
	const t_length<Sys, Y> f = focal_len<Sys, Y>(lenBefore, lenAfter);
	const Y s = Y(units::abs(units::sin(Y(0.5)*tt)));

	const t_length<Sys, Y> curv = bVert ? Y(2)*f*s : Y(2)*f/s;

	return curv;
}


// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
/**
 * @brief disc chopper burst time
 * @param r chopper radius
 * @param L chopper window length
 * @param om chopper frequency
 * @param bCounterRot single disc or two counter-rotating discs?
 * @param bSigma burst time in sigma or fwhm?
 * @return burst time
 * @see NIMA 492, pp. 97-104 (2002), doi: https://doi.org/10.1016/S0168-9002(02)01285-8
 */
template<class Sys, class Y=double>
t_time<Sys,Y> burst_time(const t_length<Sys,Y>& r,
	const t_length<Sys,Y>& L, const t_freq<Sys,Y>& om, bool bCounterRot,
	bool bSigma=1)
{
	const Y tSig = bSigma ? FWHM2SIGMA<Y> : Y(1);
	Y tScale = bCounterRot ? Y(2) : Y(1);
	return L / (r * om * tScale) * tSig;
}


/**
 * @brief disc chopper burst time
 * @see NIMA 492, pp. 97-104 (2002), doi: https://doi.org/10.1016/S0168-9002(02)01285-8
 */
template<class Sys, class Y=double>
t_length<Sys,Y> burst_time_L(const t_length<Sys,Y>& r,
	const t_time<Sys,Y>& dt, const t_freq<Sys,Y>& om, bool bCounterRot,
	bool bSigma=1)
{
	const Y tSig = bSigma ? FWHM2SIGMA<Y> : Y(1);
	Y tScale = bCounterRot ? Y(2) : Y(1);
	return dt * r * om * tScale / tSig;
}


/**
 * @brief disc chopper burst time
 * @see NIMA 492, pp. 97-104 (2002), doi: https://doi.org/10.1016/S0168-9002(02)01285-8
 */
template<class Sys, class Y=double>
t_length<Sys,Y> burst_time_r(const t_time<Sys,Y>& dt,
	const t_length<Sys,Y>& L, const t_freq<Sys,Y>& om, bool bCounterRot,
	bool bSigma=1)
{
	const Y tSig = bSigma ? FWHM2SIGMA<Y> : Y(1);
	Y tScale = bCounterRot ? Y(2) : Y(1);
	return L / (dt * om * tScale) * tSig;
}


/**
 * @brief disc chopper burst time
 * @see NIMA 492, pp. 97-104 (2002), doi: https://doi.org/10.1016/S0168-9002(02)01285-8
 */
template<class Sys, class Y=double>
t_freq<Sys,Y> burst_time_om(const t_length<Sys,Y>& r,
	const t_length<Sys,Y>& L, const t_time<Sys,Y>& dt, bool bCounterRot,
	bool bSigma=1)
{
	const Y tSig = bSigma ? FWHM2SIGMA<Y> : Y(1);
	Y tScale = bCounterRot ? Y(2) : Y(1);
	return L / (r * dt * tScale) * tSig;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------

/**
 * @brief collimation
 * @param L length of collimator
 * @param w distance between blade
 * @param bSigma calculate sigma or fwhm?
 * @return angular divergence
 * @see (Shirane 2002), Ch. 3.3
 */
template<class Sys, class Y=double>
t_angle<Sys,Y> colli_div(const t_length<Sys,Y>& L, const t_length<Sys,Y>& w, bool bSigma=1)
{
	const Y tSig = bSigma ? FWHM2SIGMA<Y> : Y(1);
	return units::atan(w/L) * tSig;
}


/**
 * @brief collimation
 * @see (Shirane 2002), Ch. 3.3
 */
template<class Sys, class Y=double>
t_length<Sys,Y> colli_div_L(const t_angle<Sys,Y>& ang, const t_length<Sys,Y>& w, bool bSigma=1)
{
	const Y tSig = bSigma ? FWHM2SIGMA<Y> : Y(1);
	return w/units::tan(ang/tSig);
}


/**
 * @brief collimation
 * @see (Shirane 2002), Ch. 3.3
 */
template<class Sys, class Y=double>
t_length<Sys,Y> colli_div_w(const t_length<Sys,Y>& L, const t_angle<Sys,Y>& ang, bool bSigma=1)
{
	const Y tSig = bSigma ? FWHM2SIGMA<Y> : Y(1);
	return units::tan(ang/tSig) * L;
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
/**
 * @brief velocity selector
 * @return selector angular frequency
 * @see https://doi.org/10.1016/0921-4526(95)00336-8
 */
template<class Sys, class Y=double>
t_freq<Sys, Y> vsel_freq(const t_length<Sys,Y>& lam,
	const t_length<Sys,Y>& len, const t_angle<Sys,Y>& twist)
{
	// https://en.wikiversity.org/wiki/De_Broglie_wavelength
	t_wavenumber<Sys,Y> k = Y(2.)*pi<Y> / lam;
	t_velocity<Sys,Y> v_n = hbar<Y>*k / m_n<Y>;

	return v_n*twist / (len * radian<Y>);
}


/**
 * @brief velocity selector
 * @see https://doi.org/10.1016/0921-4526(95)00336-8
 */
template<class Sys, class Y=double>
t_length<Sys,Y> vsel_len(const t_length<Sys,Y>& lam,
	const t_freq<Sys, Y>& om, const t_angle<Sys,Y>& twist)
{
	// https://en.wikiversity.org/wiki/De_Broglie_wavelength
	t_wavenumber<Sys,Y> k = Y(2.)*pi<Y> / lam;
	t_velocity<Sys,Y> v_n = hbar<Y>*k / m_n<Y>;

	return v_n*twist / (om * radian<Y>);
}


/**
 * @brief velocity selector
 * @see https://doi.org/10.1016/0921-4526(95)00336-8
 */
template<class Sys, class Y=double>
t_angle<Sys,Y> vsel_twist(const t_length<Sys,Y>& lam,
	const t_freq<Sys, Y>& om, const t_length<Sys,Y>& len)
{
	// https://en.wikiversity.org/wiki/De_Broglie_wavelength
	t_wavenumber<Sys,Y> k = Y(2.)*pi<Y> / lam;
	t_velocity<Sys,Y> v_n = hbar<Y>*k / m_n<Y>;

	return  (len * om * radian<Y>) / v_n;
}


/**
 * @brief velocity selector
 * @see https://doi.org/10.1016/0921-4526(95)00336-8
 */
template<class Sys, class Y=double>
t_length<Sys,Y> vsel_lam(const t_angle<Sys,Y>& twist,
	const t_freq<Sys, Y>& om, const t_length<Sys,Y>& len)
{
	t_velocity<Sys,Y> v_n = (len * om * radian<Y>) / twist;
	t_wavenumber<Sys,Y> k = m_n<Y>*v_n/hbar<Y>;
	t_length<Sys,Y> lam = Y(2.)*pi<Y> / k;

	return lam;
}

// --------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Larmor precession
//------------------------------------------------------------------------------

/**
 * gamma*B = omega
 * @see https://en.wikipedia.org/wiki/Larmor_precession
 */
template<class Sys, class Y=double>
t_freq<Sys,Y> larmor_om(const t_flux<Sys,Y>& B)
{
	return co::gamma_n * B;
}


/**
 * B = omega/gamma
 * @see https://en.wikipedia.org/wiki/Larmor_precession
 */
template<class Sys, class Y=double>
t_flux<Sys,Y> larmor_B(const t_freq<Sys,Y>& om)
{
	return om/co::gamma_n;
}


/* omega = -gamma*B
 * omega*t = -gamma*B*t
 * phi = - gamma * B * l/v
 * B = -phi*v / (gamma*l)
 * phi = -pi  =>  B = pi*v / (gamma*l)
 *
 * @see https://en.wikipedia.org/wiki/Larmor_precession
 */
template<class Sys, class Y=double>
t_flux<Sys,Y> larmor_field(const t_length<Sys,Y>& lam,
	const t_length<Sys,Y>& len,
	const t_angle<Sys,Y>& phi)
{
	t_velocity<Sys,Y> v = h<Y> / lam / co::m_n;
	t_freq<Sys,Y> om = -Y(phi/radians<Y>)*v/len;
	return om/co::gamma_n;
}

//------------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// polarisation
// ----------------------------------------------------------------------------

/**
 * polarisation density matrix
 *   (based on a proof from a lecture by P. J. Brown, 2006)
 *
 * eigenvector expansion of a state: |psi> = a_i |xi_i>
 * mean value of operator with mixed states:
 * <A> = p_i * <a_i|A|a_i>
 * <A> = tr( A * p_i * |a_i><a_i| )
 * <A> = tr( A * rho )
 * polarisation density matrix: rho = 0.5 * (1 + <P|sigma>)
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9
 * @see (Desktop Bronstein 2008), Ch. 21 (Zusatzkapitel.pdf), pp. 11-12 and p. 24
 */
template<class t_vec, class t_mat>
t_mat pol_density_mat(const t_vec& P, typename t_vec::value_type c=0.5)
requires is_vec<t_vec> && is_mat<t_mat>
{
	return (unit<t_mat>(2,2) + proj_su2<t_vec, t_mat>(P, true)) * c;
}


/**
 * Blume-Maleev equation
 * @returns scattering intensity and final polarisation vector
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9 - p. 225-226
 */
template<class t_vec, typename t_cplx = typename t_vec::value_type>
std::tuple<t_cplx, t_vec> blume_maleev(const t_vec& P_i, const t_vec& Mperp, const t_cplx& N)
requires is_vec<t_vec>
{
	const t_vec MperpConj = conj(Mperp);
	const t_cplx NConj = std::conj(N);
	constexpr t_cplx imag(0, 1);

	// ------------------------------------------------------------------------
	// intensity
	// nuclear
	t_cplx I = NConj*N;

	// nuclear-magnetic
	I += NConj*inner<t_vec>(P_i, Mperp);
	I += N*inner<t_vec>(Mperp, P_i);

	// magnetic, non-chiral
	I += inner<t_vec>(Mperp, Mperp);

	// magnetic, chiral
	I += imag * inner<t_vec>(P_i, cross<t_vec>({ MperpConj, Mperp }));
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// polarisation vector
	// nuclear
	t_vec P_f = P_i * N*NConj;

	// nuclear-magnetic
	P_f += NConj * Mperp;
	P_f += N * MperpConj;
	P_f += imag * N * cross<t_vec>({ P_i, MperpConj });
	P_f += -imag * NConj * cross<t_vec>({ P_i, Mperp });

	// magnetic, non-chiral
	P_f += Mperp * inner<t_vec>(Mperp, P_i);
	P_f += MperpConj * inner<t_vec>(P_i, Mperp);
	P_f += -P_i * inner<t_vec>(Mperp, Mperp);

	// magnetic, chiral
	P_f += imag * cross<t_vec>({ Mperp, MperpConj });
	// ------------------------------------------------------------------------

	return std::make_tuple(I, P_f/I);
}


/**
 * Blume-Maleev equation
 * calculate equation indirectly with density matrix
 *   (based on a proof from a lecture by P. J. Brown, 2006)
 *
 * V   = N*1 + <Mperp|sigma>
 * I   = tr( <V|V> rho )
 * P_f = tr( <V|sigma|V> rho ) / I
 *
 * @returns scattering intensity and final polarisation vector
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9 - p. 225-226
 */
template<class t_mat, class t_vec, typename t_cplx = typename t_vec::value_type>
std::tuple<t_cplx, t_vec> blume_maleev_indir(const t_vec& P_i, const t_vec& Mperp, const t_cplx& N)
requires is_mat<t_mat> && is_vec<t_vec>
{
	// spin-1/2
	constexpr t_cplx c = 0.5;

	// vector of pauli matrices
	const auto sigma = su2_matrices<std::vector<t_mat>>(false);

	// density matrix
	const auto density = pol_density_mat<t_vec, t_mat>(P_i, c);

	// potential
	const auto V_mag = proj_su2<t_vec, t_mat>(Mperp, true);
	const auto V_nuc = N * unit<t_mat>(2);
	const auto V = V_nuc + V_mag;
	const auto VConj = herm(V);

	// scattering intensity
	t_cplx I = trace(VConj*V * density);

	// ------------------------------------------------------------------------
	// scattered polarisation vector
	const auto m0 = (VConj * sigma[0]) * V * density;
	const auto m1 = (VConj * sigma[1]) * V * density;
	const auto m2 = (VConj * sigma[2]) * V * density;

	t_vec P_f = create<t_vec>({ trace(m0), trace(m1), trace(m2) });
	// ------------------------------------------------------------------------

	return std::make_tuple(I, P_f/I);
}


/**
 * general structure factor calculation
 * e.g. type T as vector (complex number) for magnetic (nuclear) structure factor
 * Ms_or_bs:
	- nuclear scattering lengths for nuclear neutron scattering or
	- atomic form factors for x-ray scattering
	- magnetisation (* magnetic form factor) for magnetic neutron scattering
 * Rs: atomic positions
 * Q: scattering vector G for nuclear scattering or G+k for magnetic scattering with propagation vector k
 * fs: optional magnetic form factors
 *
 * @see (Shirane 2002), p. 25, equ. 2.26 for nuclear structure factor
 * @see (Shirane 2002), p. 40, equ. 2.81 for magnetic structure factor
 * @see https://doi.org/10.1016/B978-044451050-1/50002-1
 */
template<class t_vec, class T = t_vec, template<class...> class t_cont = std::vector,
	class t_cplx = std::complex<double>>
T structure_factor(const t_cont<T>& Ms_or_bs, const t_cont<t_vec>& Rs, const t_vec& Q, const t_vec* fs=nullptr)
requires is_basic_vec<t_vec>
{
	using t_real = typename t_cplx::value_type;
	constexpr t_cplx cI{0,1};
	constexpr t_real twopi = pi<t_real> * t_real{2};
	constexpr t_real expsign = -1;

	T F{};
	if(Rs.size() == 0)
		return F;
	if constexpr(is_vec<T>)
		F = zero<T>(Rs.begin()->size());	// always 3 dims...
	else if constexpr(is_complex<T>)
		F = T(0);

	auto iterM_or_b = Ms_or_bs.begin();
	auto iterR = Rs.begin();
	typename t_vec::const_iterator iterf;
	if(fs) iterf = fs->begin();

	while(iterM_or_b != Ms_or_bs.end() && iterR != Rs.end())
	{
		// if form factors are given, use them, otherwise set to 1
		t_real f = t_real(1);
		if(fs)
		{
			auto fval = *iterf;
			if constexpr(is_complex<decltype(fval)>)
				f = fval.real();
			else
				f = fval;
		}

		// structure factor
		F += (*iterM_or_b) * f * std::exp(expsign * cI * twopi * inner<t_vec>(Q, *iterR));

		// next M or b if available (otherwise keep current)
		auto iterM_or_b_next = std::next(iterM_or_b, 1);
		if(iterM_or_b_next != Ms_or_bs.end())
			iterM_or_b = iterM_or_b_next;

		if(fs)
		{
			// next form factor if available (otherwise keep current)
			auto iterf_next = std::next(iterf, 1);
			if(iterf_next != fs->end())
				iterf = iterf_next;
		}

		// next atomic position
		std::advance(iterR, 1);
	}

	return F;
}
// ----------------------------------------------------------------------------

}
#endif
