/**
 * S(Q,w) julia interface
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2016
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#ifndef __SQWMOD_JL_H__
#define __SQWMOD_JL_H__

#include "../sqwbase.h"
#include <mutex>
#include <string>
#include <memory>


class SqwJl : public SqwBase
{
protected:
	mutable std::shared_ptr<std::mutex> m_pmtx;

	/*jl_function_t*/ void *m_pInit = nullptr;
	/*jl_function_t*/ void *m_pSqw = nullptr;
	/*jl_function_t*/ void *m_pDisp = nullptr;
	/*jl_function_t*/ void *m_pBackground = nullptr;

	// filter variables that don't start with the given prefix
	std::string m_strVarPrefix = "g_";

protected:
	std::string GetJlString(void* pVal) const;
	void PrintExceptions() const;

public:
	SqwJl() = default;
	SqwJl(const std::string& pcFile);
	virtual ~SqwJl();

	virtual std::tuple<std::vector<t_real_reso>, std::vector<t_real_reso>> disp(
		t_real_reso dh, t_real_reso dk, t_real_reso dl) const override;
	virtual t_real_reso operator()(
		t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override;
	virtual t_real_reso GetBackground(
		t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override;

	virtual std::vector<SqwBase::t_var> GetVars() const override;
	virtual void SetVars(const std::vector<SqwBase::t_var>&) override;

	virtual SqwBase* shallow_copy() const override;

	void SetVarPrefix(const char* pcFilter) { m_strVarPrefix = pcFilter; }
};

#endif
