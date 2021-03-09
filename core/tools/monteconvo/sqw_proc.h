/**
 * delegates S(Q,w) functions through processes
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2016
 * @license GPLv2
 */

#ifndef __SQW_PROC_H__
#define __SQW_PROC_H__

#include "sqw.h"
#include <mutex>
#include <memory>
#include <string>
#include <vector>
#include <unistd.h>
#include <boost/interprocess/ipc/message_queue.hpp>


enum class SqwProcStartMode
{
	// run the parent process and fork child processes
	START_PARENT_CREATE_CHILD,

	// run the parent process, forking a child process
	START_PARENT_FORK_CHILD,

	// run the child process
	START_CHILD
};


template<class t_sqw>
class SqwProc : public SqwBase
{
private:
	std::size_t m_iRefCnt = 0;


protected:
	mutable std::vector<std::shared_ptr<std::mutex>> m_pmtx;

	std::size_t m_iNumChildProcesses = 1;
	std::string m_strProcBaseName;
	std::vector<pid_t> m_pidChild;

	std::vector<std::shared_ptr<boost::interprocess::managed_shared_memory>> m_pMem;
	std::vector<std::shared_ptr<boost::interprocess::message_queue>> m_pmsgIn, m_pmsgOut;
	std::vector<void*> m_pSharedPars;


public:
	SqwProc();
	SqwProc(const char* pcCfg, SqwProcStartMode mode=SqwProcStartMode::START_PARENT_FORK_CHILD,
		const char* pcProcMemName = nullptr, const char* pcProcExecName = nullptr,
		std::size_t iNumChildProcesses = 1);
	explicit SqwProc(const std::string& strCfg);
	virtual ~SqwProc();

	virtual std::tuple<std::vector<t_real_reso>, std::vector<t_real_reso>>
		disp(t_real_reso dh, t_real_reso dk, t_real_reso dl) const override;
	virtual t_real_reso
		operator()(t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override;
	virtual bool IsOk() const override;

	virtual std::vector<SqwBase::t_var> GetVars() const override;
	virtual void SetVars(const std::vector<SqwBase::t_var>&) override;

	virtual SqwBase* shallow_copy() const override;
};

#endif
