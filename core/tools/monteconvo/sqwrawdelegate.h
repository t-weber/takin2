/**
 * interface for S(q,w) models
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 6-apr-2020
 * @license GPLv2
 */

#ifndef __MCONV_SQW_RAWDELEGATE_H__
#define __MCONV_SQW_RAWDELEGATE_H__

#include "sqwbase.h"


class SqwRawDelegate : public SqwBase
{
public:
	SqwRawDelegate() = delete;

	SqwRawDelegate(SqwBase* pDelegate) : m_pDelegate{pDelegate}
	{
	}

	virtual ~SqwRawDelegate()
	{
	}

	virtual std::tuple<std::vector<t_real_reso>, std::vector<t_real_reso>>
		disp(t_real_reso dh, t_real_reso dk, t_real_reso dl) const override
	{
		return m_pDelegate->disp(dh, dk, dl);
	}

	virtual t_real_reso operator()(t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override
	{
		return m_pDelegate->operator()(dh, dk, dl, dE);
	}

	virtual bool IsOk() const override
	{
		return m_pDelegate->IsOk();
	}


	virtual std::vector<t_var> GetVars() const override
	{
		return m_pDelegate->GetVars();
	}

	virtual const std::vector<t_var_fit>& GetFitVars() const override
	{
		return m_pDelegate->GetFitVars();
	}

	virtual void SetVars(const std::vector<t_var>& vars) override
	{
		m_pDelegate->SetVars(vars);
	}

	virtual void InitFitVars(const std::vector<t_var_fit>& vecFit) override
	{
		m_pDelegate->InitFitVars(vecFit);
	}

	virtual bool SetVarIfAvail(const std::string& strKey, const std::string& strNewVal) override
	{
		return m_pDelegate->SetVarIfAvail(strKey, strNewVal);
	}

	virtual const SqwBase& operator=(const SqwBase& sqw) override
	{
		return m_pDelegate->operator=(sqw);
	}


	virtual SqwBase* shallow_copy() const override
	{
		SqwRawDelegate *pDelegate = new SqwRawDelegate(m_pDelegate->shallow_copy());
		return pDelegate;
	}


private:
	SqwBase* m_pDelegate = nullptr;
};


#endif

