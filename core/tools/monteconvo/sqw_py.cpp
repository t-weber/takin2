/**
 * S(Q,w) python interface
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date aug-2015
 * @license GPLv2
 *
 * test:
 * g++ -DBUILD_APPLI -I. -I/usr/include/python3.7m -o tst tools/monteconvo/sqw_py.cpp tlibs/log/log.cpp tlibs/math/rand.cpp tools/monteconvo/sqwbase.cpp -lboost_system -lboost_filesystem -lboost_python37 -lpython3.7m -lrt -lpthread
 */

#include "sqw_py.h"
#include "tlibs/string/string.h"
#include "tlibs/log/log.h"
#include "tlibs/file/file.h"

#include <boost/dll/runtime_symbol_info.hpp>
#include <boost/python/stl_iterator.hpp>

using t_real = t_real_reso;

static const char* pcModIdent = "py";
static const char* pcModName = "Python Model";

#define MAX_PARAM_VAL_SIZE 128


SqwPy::SqwPy(const std::string& strFile) : m_pmtx(std::make_shared<std::mutex>())
{
	if(!tl::file_exists(strFile.c_str()))
	{
		tl::log_err("Could not find Python script file: \"", strFile, "\".");
		m_bOk = 0;
		return;
	}

	std::string strDir = tl::get_dir(strFile);
	std::string strMod = tl::get_file_noext(tl::get_file_nodir(strFile));
	const bool bSetScriptCWD = 1;

	try	// mandatory stuff
	{
		static bool bInited = 0;
		if(!bInited)
		{
			::Py_InitializeEx(0);
			if(!::Py_IsInitialized())
			{
				tl::log_err("Cannot initialise Python interpreter.");
				return;
			}

			std::string strPy = Py_GetVersion();
			tl::find_all_and_replace(strPy, std::string("\n"), std::string(", "));
			tl::log_debug("Initialised Python interpreter version ", strPy, ".");
			bInited = 1;
		}

		// set script paths
		m_sys = py::import("sys");
		py::dict sysdict = py::extract<py::dict>(m_sys.attr("__dict__"));
		py::list path = py::extract<py::list>(sysdict["path"]);
		path.append(strDir.c_str());
		path.append(".");
		std::string strSitePackages = // one directory above the takinmod_py binary
			(boost::dll::program_location().parent_path().parent_path() 
				/ "site-packages").string();
		path.append(strSitePackages.c_str());
		tl::log_debug("Using site-packages path: ", strSitePackages);

		if(bSetScriptCWD)
		{
			// set script working directory -> warning: also sets main cwd
			m_os = py::import("os");
			py::dict osdict = py::extract<py::dict>(m_os.attr("__dict__"));
			py::object pycwd = osdict["chdir"];
			if(!!pycwd)
				pycwd(strDir.c_str());
			else
				tl::log_warn("Cannot set python script working directory.");
		}

		// import takin functions
		m_mod = py::import(strMod.c_str());
		if(m_mod.is_none())
		{
			tl::log_err("Invalid Python module \"", strMod, "\".");
			m_bOk = 0;
			return;
		}

		py::dict moddict = py::extract<py::dict>(m_mod.attr("__dict__"));
		m_Sqw = moddict["TakinSqw"];
		m_bOk = !!m_Sqw;
		if(!m_bOk)
			tl::log_err("Python script has no TakinSqw function.");

		try	// optional stuff
		{
			if(moddict.has_key("TakinInit"))
			{
				m_Init = moddict["TakinInit"];
				if(!!m_Init) m_Init();
			}
			else
			{
				tl::log_warn("Python script has no TakinInit function.");
			}

			if(moddict.has_key("TakinDisp"))
				m_disp = moddict["TakinDisp"];
			else
				tl::log_warn("Python script has no TakinDisp function.");
		}
		catch(const py::error_already_set& ex) {}
	}
	catch(const py::error_already_set& ex)
	{
		PyErr_Print();
		PyErr_Clear();

		m_bOk = 0;
	}
}

SqwPy::~SqwPy()
{
	//tl::log_debug("Unloading Python interpreter.");
	//Py_FinalizeEx();
}


/**
 * E(Q)
 */
std::tuple<std::vector<t_real>, std::vector<t_real>>
	SqwPy::disp(t_real dh, t_real dk, t_real dl) const
{
	if(!m_bOk)
	{
		tl::log_err("Interpreter has not initialised, cannot query S(q,w).");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}


	std::lock_guard<std::mutex> lock(*m_pmtx);

	std::vector<t_real> vecEs, vecWs;

	try
	{
		if(!!m_disp)
		{
			py::object lst = m_disp(dh, dk, dl);
			py::stl_input_iterator<py::object> iterLst(lst);
			py::stl_input_iterator<py::object> endLst;

			if(iterLst != endLst)
			{
				py::object _vecE = *iterLst;
				py::stl_input_iterator<t_real> iterE(_vecE);
				py::stl_input_iterator<t_real> endE;

				while(iterE != endE)
					vecEs.push_back(*iterE++);
			}

			if(++iterLst != endLst)
			{
				py::object _vecW = *iterLst;
				py::stl_input_iterator<t_real> iterW(_vecW);
				py::stl_input_iterator<t_real> endW;

				while(iterW != endW)
					vecWs.push_back(*iterW++);
			}
		}
	}
	catch(const py::error_already_set& ex)
	{
		PyErr_Print();
		PyErr_Clear();
	}

	return std::make_tuple(vecEs, vecWs);
}


/**
 * S(Q,E)
 */
t_real SqwPy::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	if(!m_bOk)
	{
		tl::log_err("Interpreter has not initialised, cannot query S(q,w).");
		return t_real(0);
	}


	std::lock_guard<std::mutex> lock(*m_pmtx);
	try
	{
		return py::extract<t_real>(m_Sqw(dh, dk, dl, dE));
	}
	catch(const py::error_already_set& ex)
	{
		PyErr_Print();
		PyErr_Clear();
	}

	return 0.;
}


/**
 * Gets model variables.
 */
std::vector<SqwBase::t_var> SqwPy::GetVars() const
{
	std::vector<SqwBase::t_var> vecVars;
	if(!m_bOk)
	{
		tl::log_err("Interpreter has not initialised, cannot get variables.");
		return vecVars;
	}

	try
	{
		py::dict dict = py::extract<py::dict>(m_mod.attr("__dict__"));

		for(py::ssize_t i=0; i<py::len(dict.items()); ++i)
		{
			// name
			std::string strName = py::extract<std::string>(dict.items()[i][0]);
			if(strName.length() == 0) continue;
			if(strName[0] == '_') continue;

			// filter out non-prefixed variables
			if(m_strVarPrefix.size() && strName.substr(0,m_strVarPrefix.size()) != m_strVarPrefix)
				continue;

			// type
			std::string strType = py::extract<std::string>(dict.items()[i][1]
				.attr("__class__").attr("__name__"));
			if(strType=="module" || strType=="NoneType" || strType=="type")
				continue;
			if(strType.find("func") != std::string::npos)
				continue;

			// value
			std::string strValue = py::extract<std::string>(dict.items()[i][1].attr("__repr__")());
			if(strValue.length() > MAX_PARAM_VAL_SIZE)
			{
				//tl::log_warn("Value of variable \"", strName, "\" is too large, skipping.");
				continue;
			}

			SqwBase::t_var var;
			std::get<0>(var) = std::move(strName);
			std::get<1>(var) = std::move(strType);
			std::get<2>(var) = std::move(strValue);

			vecVars.push_back(var);
		}
	}
	catch(const py::error_already_set& ex)
	{
		PyErr_Print();
		PyErr_Clear();
	}

	return vecVars;
}


/**
 * Sets model variables.
 */
void SqwPy::SetVars(const std::vector<SqwBase::t_var>& vecVars)
{
	if(!m_bOk)
	{
		tl::log_err("Interpreter has not initialised, cannot set variables.");
		return;
	}

	try
	{
		py::dict dict = py::extract<py::dict>(m_mod.attr("__dict__"));
		std::vector<bool> vecVarsSet;
		for(std::size_t iCurVar=0; iCurVar<vecVars.size(); ++iCurVar)
			vecVarsSet.push_back(0);

		for(py::ssize_t i=0; i<py::len(dict.items()); ++i)
		{
			// variable name
			std::string strName = py::extract<std::string>(dict.items()[i][0]);
			if(strName.length() == 0) continue;
			if(strName[0] == '_') continue;

			// look for the variable name in vecVars
			bool bFound = 0;
			std::string strNewVal;
			for(std::size_t iCurVar=0; iCurVar<vecVars.size(); ++iCurVar)
			{
				const SqwBase::t_var& var = vecVars[iCurVar];
				if(std::get<0>(var) == strName)
				{
					bFound = 1;
					strNewVal = std::get<2>(var);
					vecVarsSet[iCurVar] = 1;
					break;
				}
			}
			if(!bFound)
				continue;

			// cast new value to variable type
			dict[strName] = py::eval(py::str(strNewVal), dict);
		}

		for(std::size_t iCurVar=0; iCurVar<vecVarsSet.size(); ++iCurVar)
		{
			bool bVarSet = vecVarsSet[iCurVar];
			if(!bVarSet)
			{
				tl::log_err("Could not set variable \"", std::get<0>(vecVars[iCurVar]),
					"\" as it was not found.");
			}
		}

		// TODO: check for changed parameters and if reinit is needed
		if(!!m_Init)
		{
			m_Init();
		}
	}
	catch(const py::error_already_set& ex)
	{
		PyErr_Print();
		PyErr_Clear();
	}
}


SqwBase* SqwPy::shallow_copy() const
{
	SqwPy* pSqw = new SqwPy();
	*static_cast<SqwBase*>(pSqw) = *static_cast<const SqwBase*>(this);

	pSqw->m_pmtx = this->m_pmtx;
	pSqw->m_sys = this->m_sys;
	pSqw->m_os = this->m_os;
	pSqw->m_mod = this->m_mod;
	pSqw->m_Sqw = this->m_Sqw;
	pSqw->m_Init = this->m_Init;
	pSqw->m_disp = this->m_disp;

	return pSqw;
}




// ----------------------------------------------------------------------------
// SO interface
#ifdef BUILD_PLUGIN

#include <boost/dll/alias.hpp>
#include "sqw_proc.h"
#include "sqw_proc_impl.h"
#include "libs/version.h"


std::tuple<std::string, std::string, std::string> sqw_info()
{
	//tl::log_info("In ", __func__, ".");
	return std::make_tuple(TAKIN_VER, pcModIdent, pcModName);
}


std::shared_ptr<SqwBase> sqw_construct(const std::string& strCfgFile)
{
	//tl::log_info("In ", __func__, ".");
	//return std::make_shared<SqwPy>(strCfgFile);
	return std::make_shared<SqwProc<SqwPy>>(strCfgFile);
}


SqwBase* sqw_construct_raw(const std::string& strCfgFile)
{
	return new SqwProc<SqwPy>(strCfgFile);
}


void sqw_destruct_raw(SqwBase* sqw)
{
	if(sqw) delete sqw;
}


// exports from so file
BOOST_DLL_ALIAS(sqw_info, takin_sqw_info);

// construction interface
BOOST_DLL_ALIAS(sqw_construct, takin_sqw);

// alternate raw construction interface
//BOOST_DLL_ALIAS(sqw_construct_raw, takin_sqw_new);
//BOOST_DLL_ALIAS(sqw_destruct_raw, takin_sqw_del);

#endif
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// external program
#ifdef BUILD_APPLI

#include "sqw_proc.h"
#include "sqw_proc_impl.h"
#include "libs/version.h"


int main(int argc, char** argv)
{
	if(argc <= 2)
	{
		std::cout << "#\n# This is a Takin plugin module.\n#\n";
		std::cout << "module_ident: " << pcModIdent << "\n";
		std::cout << "module_name: " << pcModName << "\n";
		std::cout << "module_type: sqw\n";
		std::cout << "required_takin_version: " << TAKIN_VER << "\n";
		std::cout.flush();
		return 0;
	}

	const char* pcCfgFile = argv[1];
	const char* pcSharedMem = argv[2];
	SqwProc<SqwPy> proc(pcCfgFile, SqwProcStartMode::START_CHILD, pcSharedMem);

	return 0;
}


#endif
// ----------------------------------------------------------------------------
