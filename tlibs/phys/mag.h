/**
 * magnetism
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 7-jul-2015 -- 2018
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_MAGDISP_H__
#define __TLIBS_MAGDISP_H__

#include <initializer_list>
#include <vector>
#include <tuple>
#include <cmath>
#include <complex>
#include <cassert>

#include "../math/linalg.h"
#include "../math/rand.h"
#include "atoms.h"
#include "nn.h"


namespace tl {
// ----------------------------------------------------------------------------

/**
 * Simple ferromagnetic dispersion
 * @param lstNeighbours list of distances to neighbour atoms and their coupling constants
 * @param vecq q position
 * @param tS spin
 * @return E(q)
 * @see e.g. (Squires 2012) p. 161
 */
template<class t_vec = ublas::vector<double>,
	typename T = typename t_vec::value_type,
	template<class...> class t_cont = std::vector>
T ferromag(const t_cont<t_vec>& vecNeighbours, const t_cont<std::complex<T>>& vecJ,
	const ublas::vector<T>& vecq, T tS)
{
	std::complex<T> J(0., 0.), J0(0., 0.);
	J = structfact<T, std::complex<T>, t_vec, t_cont>
		(vecNeighbours, vecq, vecJ, &J0).real();
	return T(2)*tS*(J0 - J).real();
}

template<class t_vec = ublas::vector<double>,
	typename T = typename t_vec::value_type,
	typename t_cont = std::initializer_list<std::tuple<t_vec, std::complex<T>>>>
T ferromag(const t_cont& lstNeighbours, const ublas::vector<T>& vecq, T tS)
{
	return ferromag(vec_from_pairvec<0,std::vector,t_cont>()(lstNeighbours),
		vec_from_pairvec<1,std::vector,t_cont>()(lstNeighbours),
		vecq, tS);
}
// ----------------------------------------------------------------------------


/**
 * Magnetic form factors
 * @see (ILL Neutron Data Booklet), sec. 2.5-1 (p. 60)
 * @see https://www.ill.eu/sites/ccsl/ffacts/
 */
template<class T=double, template<class...> class t_vec=std::initializer_list>
T j0_avg(T Q, const t_vec<T>& A, const t_vec<T>& a)
{
	assert(A.size() == a.size()+1);

	std::vector<T> vecA = A, vecB = a;
	T c = *vecA.rbegin();
	vecA.pop_back();

	// same formula as for atomic form factor, just different coefficients
	return tl::formfact<T, std::vector>(Q, vecA, vecB, c);
}

template<class T=double, template<class...> class t_vec=std::initializer_list>
T j2_avg(T Q, const t_vec<T>& A, const t_vec<T>& a)
{
	T s = Q/(T(4)*get_pi<T>());
	return j0_avg<T, t_vec>(Q, A, a) * s * s;
}

template<class T=double, template<class...> class t_vec=std::initializer_list>
T mag_formfact(T Q, T L, T S,
	const t_vec<T>& A0, const t_vec<T>& a0,
	const t_vec<T>& A2, const t_vec<T>& a2)
{
	return (L+T(2)*S) * j0_avg<T, t_vec>(Q, A0, a0) * L * j2_avg<T, t_vec>(Q, A2, a2);
}

/**
 * form factor for transition metals (d orbitals, weak LS, spin-only)
 * @see (Squires 2012), p. 138
 * @see http://www.neutron.ethz.ch/research/resources/magnetic-form-factors.html
 */
template<class T=double, template<class...> class t_vec=std::initializer_list>
std::tuple<T,T,T> mag_formfact_d(T Q, T g,
	const t_vec<T>& A0, const t_vec<T>& a0,
	const t_vec<T>& A2, const t_vec<T>& a2)
{
	T j0 = j0_avg<T, t_vec>(Q, A0, a0);
	T j2 = j2_avg<T, t_vec>(Q, A2, a2);

	return std::tuple<T,T,T>(j0, j2, j0 + (T(1)-T(2)/g)*j2);
}

/**
 * form factor for rare earths (f orbitals, strong LS, jj)
 * @see (Squires 2012), p. 139
 * @see http://www.neutron.ethz.ch/research/resources/magnetic-form-factors.html
 */
template<class T=double, template<class...> class t_vec=std::initializer_list>
std::tuple<T,T,T> mag_formfact_f(T Q, T L, T S, T J,
	const t_vec<T>& A0, const t_vec<T>& a0,
	const t_vec<T>& A2, const t_vec<T>& a2)
{
	T j0 = j0_avg<T, t_vec>(Q, A0, a0);
	T j2 = j2_avg<T, t_vec>(Q, A2, a2);

	T gL = T(0.5) + (L*(L+T(1)) - S*(S+T(1))) / (T(2)*J*(J+T(1)));
	T gS = T(1) + (S*(S+T(1)) - L*(L+T(1))) / (J * (J+T(1)));

	return std::tuple<T,T,T>(j0, j2, (gS*j0 + gL*(j0+j2)) / (gL + gS));
}
// ----------------------------------------------------------------------------


/**
 * effective magnetic scattering length p in fm
 * @param tM: magnetic moment (in units of muB)
 * @see e.g.: (Shirane 2002), p. 4
 */
template<typename T = double>
T mag_scatlen_eff(T tM)
{
	const T cFact = tl::get_r_e<T>() *
		(-tl::get_mu_n<T>()/tl::get_mu_N<T>()) * T(0.5)
		/ (tl::get_one_femtometer<T>());

	return cFact * tM;
}


/**
 * spin S perpendicular to scattering vector Q
 * @see (Shirane 2002), p. 37, equ. 2.63
 */
template<class t_vec = ublas::vector<double>>
t_vec get_S_perp_Q(const t_vec& S, const t_vec& Q)
{
	t_vec Qnorm = Q / veclen(Q);
	// subtract parts of S not perpendicular to Q
	return S - mult<t_vec, t_vec>(Qnorm, S)*Qnorm;
}


/**
 * calculates the magnetic structure factor Fm
 * @param lstAtoms list of atom positions
 * @param lstAtoms list of spins
 * @param vecQ scattering vector
 * @param lstf G-dependent magnetic form factor
 * @param lstg g factors, 2 if none given
 * @param pF0 optional total form factor.
 * @param dVuc optionally normalise by the unit cell volume
 * @return structure factor
 * @see (Shirane 2002), p. 40, equ. 2.81
 */
template<typename T = double, typename t_ff = std::complex<T>,
	template<class...> class t_vec = ublas::vector,
	template<class ...> class t_cont = std::initializer_list>
t_vec<std::complex<T>> structfact_mag(const t_cont<t_vec<T>>& lstAtoms,
	const t_cont<t_vec<T>>& lstS, const t_vec<T>& vecQ,
	const t_cont<t_ff>& lstf = t_cont<t_ff>(),
	const t_cont<T>& lstg = t_cont<T>(),
	t_ff *pF0 = nullptr, T dVuc = T(-1))
{
	constexpr std::complex<T> i(0., 1.);
	t_vec<std::complex<T>> Fm = make_vec<t_vec<std::complex<T>>>({{0., 0.},{0., 0.},{0., 0.}});

	using t_iter_atoms = typename t_cont<t_vec<T>>::const_iterator;
	using t_iter_ffact = typename t_cont<t_ff>::const_iterator;
	using t_iter_spin = typename t_cont<t_vec<T>>::const_iterator;
	using t_iter_g = typename t_cont<T>::const_iterator;

	t_iter_atoms iterAtom = lstAtoms.begin();
	t_iter_ffact iterFFact = lstf.begin();
	t_iter_spin iterSpin = lstS.begin();
	t_iter_g iterg = lstg.begin();

	/*const T cFact = tl::get_r_e<T>() * tl::get_gamma_n<T>() * T(0.5) *
		tl::get_one_second<T>() * tl::get_one_tesla<T>()
		/ (tl::get_one_meter<T>() * T(1e-15));*/
	const T cFact = tl::get_r_e<T>() *
		(-tl::get_mu_n<T>()/tl::get_mu_N<T>()) * T(0.5)
		/ (tl::get_one_meter<T>() * T(1e-15));	// in fm

	if(pF0) *pF0 = t_ff(0);

	for(; iterAtom != lstAtoms.end(); ++iterAtom)
	{
		// only use form factors when available
		t_ff tFF = T(1);
		if(iterFFact != lstf.end())
			tFF = *iterFFact;

		// use g=2 if none given
		T g = T(2);
		if(iterg != lstg.end())
			g = *iterg;

		const t_vec<T>& vecS = *iterSpin;
		t_vec<std::complex<T>> vecSperp = get_S_perp_Q<t_vec<T>>(vecS, vecQ);

		Fm += cFact * g * tFF * vecSperp *
			std::exp(i * (mult<t_vec<T>, t_vec<T>>(vecQ, *iterAtom)));
		if(pF0) *pF0 += tFF;

		// if there is only one form factor in the list, use it for all positions
		if(iterFFact!=lstf.end() && std::next(iterFFact)!=lstf.end())
			++iterFFact;
		// if there is only one g factor in the list, use it for all positions
		if(iterg!=lstg.end() && std::next(iterg)!=lstg.end())
			++iterg;
		// if there is only one spin in the list, use it for all positions
		if(iterSpin!=lstS.end() && std::next(iterSpin)!=lstS.end())
			++iterSpin;
	}

	// if unit cell volume is given, normalise by it
	if(dVuc > T(0)) Fm /= dVuc;
	return Fm;
}




// ----------------------------------------------------------------------------
/**
 * metropolis algorithm
 * @see (Scherer 2010), p. 104
 * @see (Schroeder 2000), p. 346 ff
 */
template<class t_real, std::size_t DIM,
	template<class, std::size_t, class...> class t_arr_1d = boost::array,
	template<class, std::size_t, class...> class t_arr_nd = boost::multi_array>
void metrop(
	const t_arr_1d<typename t_arr_nd<bool,DIM>::index, DIM>& arrDims,
	std::size_t iNumIters, t_real dJ, t_real dk, t_real dT,
	t_arr_nd<bool, DIM>& arrSpins,
	t_real *pEtot = nullptr, t_real *pS = nullptr)
{
	using T = bool;
	using t_arr = t_arr_nd<T, DIM>;
	using t_idx = typename t_arr::index;
	using t_dim = t_arr_1d<t_idx, DIM>;

	const t_real dBeta = t_real(1)/(dk*dT);

	// get next neighbours
	auto getNN = [](const t_dim& dim, const t_dim& idx) -> std::vector<t_dim>
	{
		std::vector<t_dim> vecNN;
		for(std::size_t iDim=0; iDim<DIM; ++iDim)
		{
			if(idx[iDim] > 0)
			{
				t_dim idxNew = idx;
				--idxNew[iDim];
				vecNN.emplace_back(std::move(idxNew));
			}
			if(idx[iDim] < dim[iDim]-1)
			{
				t_dim idxNew = idx;
				++idxNew[iDim];
				vecNN.emplace_back(std::move(idxNew));
			}
		}
		return vecNN;
	};

	// calculate energy
	auto calcE = [dJ](const t_arr& arr,
		const t_dim& idxSpin, const std::vector<t_dim>& vecNN,
		bool bFlip) -> t_real
	{
		t_real dE = t_real(0);
		bool bSpin = arr(idxSpin);
		if(bFlip) bSpin = !bSpin;
		t_real dSpin = (bSpin ? t_real(1) : t_real(-1));

		for(const t_dim& idxNN : vecNN)
		{
			t_real dSpinNN = (arr(idxNN) ? t_real(1) : t_real(-1));
			dE += -dJ*dSpin*dSpinNN;
		}
		return dE;
	};


	// TODO: search radius
	t_dim dimMin; dimMin.fill(0);
	t_dim dimMax = arrDims;

	//t_arr arrRand = rand_array<T, DIM, t_arr_1d, t_arr_nd>(arrDims);
	//t_arr* parrSpins = &arrRand;
	t_arr *parrSpins = &arrSpins;

	for(std::size_t iIter=0; iIter<iNumIters; ++iIter)
	{
		t_dim idx = rand_idx<t_idx, DIM, t_arr_1d>(dimMin, dimMax);
		std::vector<t_dim> vecNN = getNN(arrDims, idx);

		t_real dENoFlip = calcE(*parrSpins, idx, vecNN, 0);
		t_real dEFlip = calcE(*parrSpins, idx, vecNN, 1);

		if(dEFlip < dENoFlip)
		{
			(*parrSpins)(idx) = !(*parrSpins)(idx);
		}
		else
		{
			t_real dEDiff = dEFlip - dENoFlip;
			t_real dProb = std::exp(-dBeta * dEDiff);

			if(rand_prob<t_real>(dProb))
				(*parrSpins)(idx) = !(*parrSpins)(idx);
		}
	}

	if(pEtot)	// calculate total energy
	{
		*pEtot = t_real(0);

		t_dim idxCur; idxCur.fill(0);
		while(1)
		{
			std::vector<t_dim> vecNN = getNN(arrDims, idxCur);
			t_real dENoFlip = calcE(*parrSpins, idxCur, vecNN, 0);
			//t_real dEFlip = calcE(*parrSpins, idxCur, vecNN, 1);

			//dZ += dBoltzNoFlip + dBoltzFlip;	// partition function
			*pEtot += dENoFlip;

			++idxCur[0];
			for(std::size_t iDim=0; iDim<DIM-1; ++iDim)
			{
				if(idxCur[iDim] >= arrDims[iDim])
				{
					idxCur[iDim] = 0;
					++idxCur[iDim+1];
				}
			}
			if(idxCur[DIM-1] >= arrDims[DIM-1])
				break;
		}

		// normalise
		*pEtot /= std::pow(2., t_real(DIM));		// NN
		for(std::size_t iDim=0; iDim<DIM; ++iDim)	// N
			*pEtot /= t_real(arrDims[iDim]);
	}

	if(pS)	// calculate mean magnetisation
	{
		*pS = t_real(0);

		const T* pDat = parrSpins->data();
		std::size_t iNumElems = parrSpins->num_elements();

		for(std::size_t iElem=0; iElem<iNumElems; ++iElem)
			*pS += (*(pDat+iElem) ? t_real(1) : t_real(-1));

		// normalise
		*pS /= std::pow(2., t_real(DIM));		// NN
		for(std::size_t iDim=0; iDim<DIM; ++iDim)	// N
			*pS /= t_real(arrDims[iDim]);
	}
}



// ----------------------------------------------------------------------------
// special paramagnetic functions
// ----------------------------------------------------------------------------

/**
 * Langevin function for mean cosine
 * @see https://en.wikipedia.org/wiki/Brillouin_and_Langevin_functions
 */
template<class T=double>
T langevin(T x)
{
	return tl::coth(x) - T(1)/x;
}


/**
 * Brillouin function ~ magnetisation
 * @param x = g muB J B / (kB T)
 * @see https://en.wikipedia.org/wiki/Brillouin_and_Langevin_functions
 */
template<class T=double>
T brillouin(T J, T x)
{
	T Jfact = T(1)+T(1)/(T(2)*J);
	T Jfact2 = T(1)/(T(2)*J);

	return Jfact * tl::coth(Jfact*x)
		- Jfact2 * tl::coth(Jfact2*x);
}

// ----------------------------------------------------------------------------


}
#endif
