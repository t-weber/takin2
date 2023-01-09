/**
 * dummy S(q,w) model
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 8-apr-2020
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

#ifndef __MCONV_SQW_DUMMY_H__
#define __MCONV_SQW_DUMMY_H__

#include "sqwbase.h"


class SqwNull : public SqwBase
{
public:
	SqwNull() {}
	SqwNull(const char*) {}
	~SqwNull() {}

	virtual std::tuple<std::vector<t_real_reso>, std::vector<t_real_reso>>
		disp(t_real_reso dh, t_real_reso dk, t_real_reso dl) const override
	{
		return std::make_tuple(std::vector<t_real_reso>{}, std::vector<t_real_reso>{});
	}

	virtual t_real_reso operator()(t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override
	{
		return t_real_reso{0};
	}

	virtual bool IsOk() const override
	{
		return 1;
	}

	virtual std::vector<t_var> GetVars() const override
	{
		return std::vector<t_var>{};
	}

	virtual const std::vector<t_var_fit>& GetFitVars() const override
	{
		static std::vector<t_var_fit> vec;
		return vec;
	}

	virtual void SetVars(const std::vector<t_var>& vars) override
	{
	}

	virtual void InitFitVars(const std::vector<t_var_fit>& vecFit) override
	{
	}

	virtual bool SetVarIfAvail(const std::string& strKey, const std::string& strNewVal) override
	{
		return 1;
	}

	virtual const SqwBase& operator=(const SqwBase& sqw) override
	{
		return *this;
	}

	virtual SqwBase* shallow_copy() const override
	{
		return new SqwNull();
	}

};


#endif

