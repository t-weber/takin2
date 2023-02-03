/**
 * Powder peaks
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date apr-2015
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
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

#ifndef __POWDER_H__
#define __POWDER_H__

#include <tuple>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <initializer_list>
#include "lattice.h"
#include "../string/string.h"
#include "../helper/hash.h"
#include "../helper/array.h"
#include "../fit/minuit.h"
#include "../phys/neutrons.h"


namespace tl {

template<class t_int=int, class t_real=double>
class Powder
{
	public:
		// hkl G F
		typedef std::tuple<t_int, t_int, t_int, t_real, t_real> t_peak;
		typedef std::string t_str;

	private:
		static t_str to_str(t_real t)
		{
			static const int s_iPrec = 8;
			return tl::var_to_str<t_real, t_str>(t, s_iPrec);
		}


		static bool is_eq(t_real t0, t_real t1)
		{
			t_str str0 = to_str(t0);
			t_str str1 = to_str(t1);
			return str0 == str1;
		}


		static size_t hash_peak(const t_peak& peak)
		{
			return tl::hash_ordered<std::initializer_list<t_int>>
			({
				std::get<0>(peak),
				std::get<1>(peak),
				std::get<2>(peak)
			});
		}


		static size_t hash_peak_unique(const t_peak& peak)
		{
			t_str strG = to_str(std::get<3>(peak));
			t_str strF = to_str(std::get<4>(peak));

			return tl::hash_ordered<std::initializer_list<std::string>>({strG, strF});
		}


		static bool equ_peak(const t_peak& peak0, const t_peak& peak1)
		{
			return std::get<0>(peak0) == std::get<0>(peak1) &&
				std::get<1>(peak0) == std::get<1>(peak1) &&
				std::get<2>(peak0) == std::get<2>(peak1);
		}


		static bool equ_peak_unique(const t_peak& peak0, const t_peak& peak1)
		{
			bool bGEq = is_eq(std::get<3>(peak0), std::get<3>(peak1));
			bool bFEq = is_eq(std::get<4>(peak0), std::get<4>(peak1));

			return bGEq && bFEq;
		}

	public:
		// all peaks
		typedef std::unordered_set<t_peak, decltype(&hash_peak), decltype(&equ_peak)> t_peaks;
		// peaks unique in G and F
		typedef std::unordered_set<t_peak, decltype(&hash_peak_unique), decltype(&equ_peak_unique)> t_peaks_unique;

	protected:
		t_peaks m_peaks;					// hashes & compares hkl
		t_peaks_unique m_peaks_unique;		// hashes & compares F & G

		// associated reciprocal lattice
		const Lattice<t_real> *m_pLatticeRecip = nullptr;

	public:
		Powder()
			: m_peaks(10, &hash_peak, &equ_peak),
			  m_peaks_unique(10, &hash_peak_unique, &equ_peak_unique)
		{}


		ublas::vector<t_real> GetRecipLatticePos(t_real dh, t_real dk, t_real dl) const
		{
			if(m_pLatticeRecip)
				return m_pLatticeRecip->GetPos(dh, dk, dl);

			return ublas::vector<t_real>();
		}


		t_real GetG(t_real dh, t_real dk, t_real dl) const
		{
			ublas::vector<t_real> vecG = GetRecipLatticePos(dh, dk, dl);
			if(vecG.size())
				return veclen(vecG);
			return 0.;
		}


		void AddPeak(t_int h, t_int k, t_int l, t_real F=0.)
		{
			t_peak peak(h,k,l, GetG(h,k,l), F);

			m_peaks.insert(peak);
			m_peaks_unique.insert(peak);
		}


		void SetRecipLattice(const Lattice<t_real>* pLatt)
		{
			m_pLatticeRecip = pLatt;
		}

		const t_peaks& GetPeaks() const { return m_peaks; }
		const t_peaks& GetUniquePeaks() const { return m_peaks_unique; }


		bool HasPeak(int h, int k, int l) const
		{
			t_peak peak(h,k,l);
			return m_peaks.find(peak) != m_peaks.end();
		}


		bool HasUniquePeak(int h, int k, int l, t_real F) const
		{
			const t_real G = GetG(h,k,l);

			for(const t_peak& pk : m_peaks_unique)
			{
				if(is_eq(G, std::get<3>(pk)) && is_eq(F, std::get<4>(pk)))
					return 1;
			}
			return 0;
		}


		std::size_t GetMultiplicity(t_int h, t_int k, t_int l) const
		{
			t_real G = GetG(h,k,l);

			std::size_t iMult = 0;
			for(const t_peak& pk : m_peaks)
			{
				if(is_eq(G, std::get<3>(pk)))
					++iMult;
			}
			return iMult;
		}


		std::size_t GetMultiplicity(t_int h, t_int k, t_int l, t_real F) const
		{
			t_real G = GetG(h,k,l);

			std::size_t iMult = 0;
			for(const t_peak& pk : m_peaks)
			{
				if(is_eq(G, std::get<3>(pk)) && is_eq(F, std::get<4>(pk)))
					++iMult;
			}
			return iMult;
		}


		/**
		 * get peaks only unique in G (not F)
		 */
		std::vector<t_peak> GetUniquePeaksSumF() const
		{
			std::vector<t_peak> vecPeaks;

			for(const t_peak& peak : m_peaks)
			{
				typename decltype(vecPeaks)::iterator iter = std::find_if(
					vecPeaks.begin(), vecPeaks.end(), [&peak](const t_peak& pk) -> bool
				{
					t_real curG = std::get<3>(peak);
					t_real pkG = std::get<3>(pk);

					return is_eq(curG, pkG);
				});

				// not already in vector?
				if(iter == vecPeaks.end())
					vecPeaks.push_back(peak);
				else
					std::get<4>(*iter) += std::get<4>(peak);	// add F
			}

			return vecPeaks;
		}


		void clear()
		{
			m_peaks.clear();
			m_peaks_unique.clear();
			m_pLatticeRecip = nullptr;
		}
};



/**
 * corrects mono & sample axes using known powder lines
 * @desc see (Shirane 2002), p. 87
 */
template<class t_real = double>
bool powder_align(t_real _d_mono, const std::vector<t_real>& vecGs,
	const std::vector<t_real>& vecTTs, const std::vector<t_real>& vecTTErrs,
	std::vector<t_real>& vecRes, std::vector<t_real>& vecResErrs)
{
	// fit function
	auto fktBragg = [](t_real _G, t_real _k, t_real _dtt) -> t_real
	{
		t_wavenumber_si<t_real> G = _G / get_one_angstrom<t_real>();
		t_wavenumber_si<t_real> k = _k / get_one_angstrom<t_real>();
		t_angle_si<t_real> twotheta;

		try
		{
			twotheta = bragg_recip_twotheta(G, k, t_real(1));
		}
		catch(const std::exception& ex)
		{
			tl::log_err(ex.what());
			return t_real(std::numeric_limits<t_real>::max());
		}

		return t_real(twotheta/get_one_radian<t_real>()) + _dtt;
	};


	// nominal ki given as fitting hint
	t_real dNominalKi = 0.;
	if(vecRes.size() >= 1) dNominalKi = vecRes[0];


	// conversion
	std::vector<t_real_min> _vecRes =
		container_cast<t_real_min, t_real, std::vector>()(vecRes);
	std::vector<t_real_min> _vecResErrs =
		container_cast<t_real_min, t_real, std::vector>()(vecResErrs);

	std::vector<std::string> vecParams = { "k", "dtt" };
	bool bFitOk = fit<3>(fktBragg,
		container_cast<t_real_min, t_real, std::vector>()(vecGs), 	// x
		container_cast<t_real_min, t_real, std::vector>()(vecTTs),	// y
		container_cast<t_real_min, t_real, std::vector>()(vecTTErrs),	// yerr
		vecParams, _vecRes, _vecResErrs);

	// back-conversion
	vecRes = container_cast<t_real, t_real_min, std::vector>()(_vecRes);
	vecResErrs = container_cast<t_real, t_real_min, std::vector>()(_vecResErrs);


	// calc. mono angle
	if(vecRes.size() >= 3)
	{
		t_wavenumber_si<t_real> ki = vecRes[0] / get_one_angstrom<t_real>();
		t_length_si<t_real> d_mono = _d_mono * get_one_angstrom<t_real>();
		t_angle_si<t_real> twotheta_m = bragg_recip_twotheta(d2G(d_mono), ki, t_real(1));
		vecRes[2] = t_real(twotheta_m / get_one_radian<t_real>());
		vecResErrs[2] = t_real(-1);	// TODO: propagate error
	}

	// calc. mono angle delta
	if(vecRes.size() >= 4)
	{
		t_wavenumber_si<t_real> ki = dNominalKi / get_one_angstrom<t_real>();
		t_length_si<t_real> d_mono = _d_mono * get_one_angstrom<t_real>();
		t_angle_si<t_real> twotheta_m = bragg_recip_twotheta(d2G(d_mono), ki, t_real(1));
		t_real dNominalAngle = t_real(twotheta_m / get_one_radian<t_real>());
		vecRes[3] = vecRes[2] - dNominalAngle;
		vecResErrs[3] = t_real(-1);	// TODO: propagate error
	}

	return bFitOk;
}


}

#include <ostream>

template<class t_int=int, class t_real=double>
std::ostream& operator<<(std::ostream& ostr, const tl::Powder<t_int,t_real>& powder)
{
	const typename tl::Powder<t_int,t_real>::t_peaks& peaks = powder.GetPeaks();
	const typename tl::Powder<t_int,t_real>::t_peaks_unique& peaks_unique = powder.GetUniquePeaks();

	ostr << "Peaks:\n";
	for(const typename tl::Powder<t_int,t_real>::t_peak& pk : peaks)
	{
		t_int h = std::get<0>(pk);
		t_int k = std::get<1>(pk);
		t_int l = std::get<2>(pk);

		ostr << "\t(" << h << k << l << ")\n";
	}

	ostr << "Unique Peaks:\n";
	for(const typename tl::Powder<t_int,t_real>::t_peak& pk : peaks_unique)
	{
		t_int h = std::get<0>(pk);
		t_int k = std::get<1>(pk);
		t_int l = std::get<2>(pk);

		ostr << "\t(" << h << k << l << ")";
		ostr << ", multiplicity: " << powder.GetMultiplicity(h,k,l) << "\n";
	}

	return ostr;
}

#endif
