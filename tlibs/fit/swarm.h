/**
 * Swarm-fitting algorithms
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2017
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

#ifndef __SWARM_H__
#define __SWARM_H__

#include <vector>
#include <functional>

#include "funcmod.h"
#include "../math/linalg.h"
#include "../math/math.h"
#include "../math/rand.h"
#include "../log/log.h"


namespace tl {


/**
 * a single swarm particle
 */
template<class t_real, template<class...> class t_vec>
struct Raven
{
	t_vec<t_real> vecPos, vecBestPos;
	t_vec<t_real> vecVel;

	void TimeStep(t_real dDelta)
	{
		vecPos += dDelta*vecVel;
	}
};



/**
 * swarm minimisation
 * algorithm: https://en.wikipedia.org/wiki/Particle_swarm_optimization
 */
template<class t_real, template<class...> class t_vec>
class Unkindness
{
protected:
	std::vector<Raven<t_real, t_vec>> m_vecRavens;
	t_vec<t_real> m_vecBestPos;
	bool m_bBestPos = 0;

	// function to fit
	FitterFuncModel<t_real> *m_pFunc = nullptr;

	// measured 1
	std::size_t m_iSizeDat = 0;
	const t_real *m_pDatX = nullptr;
	const t_real *m_pDatY = nullptr;
	const t_real *m_pDatYErr = nullptr;

	t_real m_dVelScale = 1.;
	t_real m_dPartScale = 1.2;
	t_real m_dSwarmScale = 1.4;
	t_real m_dTimeDelta = 1;
	t_real m_dEps = get_epsilon<t_real>();

	std::size_t m_iMaxCalls = 500;			// max. func. calls
	std::size_t m_iMaxBestPosIters = 5;		// max. iterations without improvement


public:
	void SetData(std::size_t iSize,
		const t_real *pX, const t_real *pY, const t_real *pYErr)
	{
		m_iSizeDat = iSize;

		m_pDatX = pX;
		m_pDatY = pY;
		m_pDatYErr = pYErr;
	}


	// ------------------------------------------------------------------------
	void SetMaxCalls(std::size_t iMaxCalls) { m_iMaxCalls = iMaxCalls; }
	void SetMaxBestIters(std::size_t iMaxIters) { m_iMaxBestPosIters = iMaxIters; }

	void SetVelScale(t_real dSc) { m_dVelScale = dSc; }
	void SetPartScale(t_real dSc) { m_dPartScale = dSc; }
	void SetSwarmScale(t_real dSc) { m_dSwarmScale = dSc; }
	void SetTimeDelta(t_real dSc) { m_dTimeDelta = dSc; }

	void SetEpsilon(t_real dEps) { m_dEps = dEps; }

	const t_vec<t_real>& GetBestPos() const { return m_vecBestPos; }
	bool IsBestPosValid() const { return m_bBestPos; }
	// ------------------------------------------------------------------------


	t_real EvalFunc(const t_vec<t_real>& vecParams)
	{
		if(!m_pFunc) return t_real(0);
		m_pFunc->SetParams(convert_vec_full<t_real, t_real, t_vec, std::vector>(vecParams));

		// wrapper for m_pFunc
		auto func = [this](t_real x) -> t_real
		{
			t_real y = (*m_pFunc)(x);
			//log_debug("func(", x, ") = ", y);
			return y;
		};

		t_real dchi2 = chi2<t_real, decltype(func), decltype(m_pDatX)>(
			func, m_iSizeDat, m_pDatX, m_pDatY, m_pDatYErr);
		return dchi2;
	}


	void Init(std::size_t iNumRavens,
		FitterFuncModel<t_real>* pFunc,
		const t_vec<t_real>& vecMin, const t_vec<t_real>& vecMax)
	{
		m_pFunc = pFunc;

		m_vecRavens.clear();
		m_vecRavens.reserve(iNumRavens);

		// random initial best position
		m_vecBestPos = convert_vec_full<t_real, t_real, std::vector, t_vec>(
			rand_minmax_nd<t_real, std::vector>(
			convert_vec_full<t_real, t_real, t_vec, std::vector>(vecMin),
			convert_vec_full<t_real, t_real, t_vec, std::vector>(vecMax)));

		t_vec<t_real> vecVelMin = vecMin-vecMax;
		t_vec<t_real> vecVelMax = vecMax-vecMin;

		for(std::size_t i=0; i<iNumRavens; ++i)
		{
			Raven<t_real, t_vec> raven;

			// random positions
			raven.vecPos = raven.vecBestPos =
				convert_vec_full<t_real, t_real, std::vector, t_vec>(
					rand_minmax_nd<t_real, std::vector>(
					convert_vec_full<t_real, t_real, t_vec, std::vector>(vecMin),
					convert_vec_full<t_real, t_real, t_vec, std::vector>(vecMax)));

			// random velocities
			raven.vecVel =
				convert_vec_full<t_real, t_real, std::vector, t_vec>(
					rand_minmax_nd<t_real, std::vector>(
					convert_vec_full<t_real, t_real, t_vec, std::vector>(vecVelMin),
					convert_vec_full<t_real, t_real, t_vec, std::vector>(vecVelMax)));

			if(EvalFunc(raven.vecPos) < EvalFunc(m_vecBestPos))
				m_vecBestPos = raven.vecPos;

			m_vecRavens.emplace_back(std::move(raven));
		}
	}


	void Run()
	{
		if(!m_vecRavens.size()) return;
		const std::size_t iDim = m_vecRavens[0].vecPos.size();

		m_bBestPos = 1;
		std::size_t iFunc = 0;
		t_vec<t_real> vecLastBestPos;
		std::size_t iLastBestPos = 0;

		while(1)
		{
			t_vec<t_real> vec0 = fill_vector<t_vec<t_real>>(iDim, t_real(0));
			t_vec<t_real> vec1 = fill_vector<t_vec<t_real>>(iDim, t_real(1));

			for(Raven<t_real, t_vec>& raven : m_vecRavens)
			{
				// random vectors
				t_vec<t_real> vec01_part =
					convert_vec_full<t_real, t_real, std::vector, t_vec>(
						rand_minmax_nd<t_real, std::vector>(
						convert_vec_full<t_real, t_real, t_vec, std::vector>(vec0),
						convert_vec_full<t_real, t_real, t_vec, std::vector>(vec1)));

				t_vec<t_real> vec01_swarm =
					convert_vec_full<t_real, t_real, std::vector, t_vec>(
						rand_minmax_nd<t_real, std::vector>(
						convert_vec_full<t_real, t_real, t_vec, std::vector>(vec0),
						convert_vec_full<t_real, t_real, t_vec, std::vector>(vec1)));

				// new velocity
				raven.vecVel = m_dVelScale*raven.vecVel
					+ m_dPartScale*ublas::element_prod(vec01_part, raven.vecBestPos-raven.vecPos)
					+ m_dSwarmScale*ublas::element_prod(vec01_swarm, m_vecBestPos-raven.vecPos);

				raven.TimeStep(m_dTimeDelta);

				// new best minimum positions
				if(EvalFunc(raven.vecPos) < EvalFunc(raven.vecBestPos))
				{
					raven.vecBestPos = raven.vecPos;
					if(EvalFunc(raven.vecPos) < EvalFunc(m_vecBestPos))
						m_vecBestPos = raven.vecPos;
				}
			}

			// no new best position since a few iterations?
			if(vec_equal<t_vec<t_real>>(m_vecBestPos, vecLastBestPos, m_dEps))
				++iLastBestPos;
			else
				iLastBestPos = 0;
			if(iLastBestPos >= m_iMaxBestPosIters)
				return;

			// max. function calls reached
			if(++iFunc >= m_iMaxCalls)
			{
				tl::log_warn("Maximum number of function calls reached.");
				m_bBestPos = 0;
				return;
			}

			vecLastBestPos = m_vecBestPos;
		}
	}
};



// -----------------------------------------------------------------------------
template<typename t_real, std::size_t iNumArgs, typename t_func>
bool swarmfit(t_func&& func,
	const std::vector<t_real>& vecX,
	const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr,

	const std::vector<std::string>& vecParamNames,	// size: iNumArgs-1
	std::vector<t_real>& vecVals,
	std::vector<t_real>& vecErrs,

	bool bDebug = 1)
{
	if(!vecX.size() || !vecY.size() || !vecYErr.size())
	{
		tl::log_err("No data given to fitter.");
		return false;
	}

	std::size_t iParamSize = iNumArgs-1;
	if(iParamSize != vecVals.size())
	{
		tl::log_err("Parameter size mismatch.");
		return false;
	}

	std::size_t iDatSize = std::min(vecY.size(), vecYErr.size());
	std::size_t iNumRavens = 256*iParamSize;
	std::size_t iMaxCalls = 128;
	std::size_t iMaxBestIters = 4;

	ublas::vector<t_real> vecMin(iParamSize), vecMax(iParamSize);
	for(std::size_t iY=0; iY<iParamSize; ++iY)
	{
		vecMin[iY] = vecVals[iY] - vecErrs[iY];
		vecMax[iY] = vecVals[iY] + vecErrs[iY];
	}

	FitterLamFuncModel<t_real, iNumArgs, t_func> mod(func);
	tl::Unkindness<t_real, ublas::vector> unk;
	unk.SetData(iDatSize, vecX.data(), vecY.data(), vecYErr.data());
	unk.SetMaxCalls(iMaxCalls);
	unk.SetMaxBestIters(iMaxBestIters);
	unk.Init(iNumRavens, &mod, vecMin, vecMax);
	unk.Run();

	const auto& vecBest = unk.GetBestPos();

	// check if results are within error ranges
	bool bInRange = 1;
	for(std::size_t iParam=0; iParam<vecBest.size(); ++iParam)
	{
		if(!is_in_range<t_real>(vecBest[iParam], vecVals[iParam], vecErrs[iParam]))
		{
			bInRange = 0;
			break;
		}
	}

	// copy results
	for(std::size_t iParam=0; iParam<vecBest.size(); ++iParam)
	{
		vecVals[iParam] = vecBest[iParam];
		//vecErrs[iParam] = 0.;
	}

	if(bDebug)
	{
		std::ostringstream ostrRes;
		for(std::size_t iParam=0; iParam<vecBest.size(); ++iParam)
			ostrRes << vecParamNames[iParam] << " = " << vecVals[iParam] << ", ";

		tl::log_debug("Swarm fit: valid = ", unk.IsBestPosValid(),
			", in_range = ", bInRange,
			", result: ", ostrRes.str());
	}

	return unk.IsBestPosValid() /*&& bInRange*/;
}
// -----------------------------------------------------------------------------

}

#endif
