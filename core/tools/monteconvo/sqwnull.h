/**
 * dummy S(q,w) model
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 8-apr-2020
 * @license GPLv2
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

