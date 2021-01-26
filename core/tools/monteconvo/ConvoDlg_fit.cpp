/**
 * monte carlo convolution tool -> convolution fitting
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015, 2016
 * @license GPLv2
 */

#include "ConvoDlg.h"

#include <QMutex>
#include <QMessageBox>

#include <vector>
#include <string>
#include <exception>

#include <boost/optional.hpp>

#include "tlibs/string/string.h"


using t_real = t_real_reso;
using t_real_min = double;


#include <Minuit2/FCNBase.h>
#include <Minuit2/MnFcn.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnPrint.h>

namespace minuit = ROOT::Minuit2;


class StopRequestedEx : public std::runtime_error
{
public:
	StopRequestedEx(const char* msg) : runtime_error{msg} {}
};



/**
  * interface between ConvoDlg and Minuit
  */
class MinuitFunc : public minuit::FCNBase
{
public:
	MinuitFunc(ConvoDlg* convodlg, const ConvoDlg::t_sqwparams* sqwparams)
		: m_convodlg{convodlg}, m_sqwparams{sqwparams}, m_seed{tl::get_rand_seed()}
	{}

	MinuitFunc(const MinuitFunc&) = delete;
	const MinuitFunc& operator=(const MinuitFunc&) = delete;

	virtual ~MinuitFunc() = default;


	virtual t_real_min operator()(const std::vector<t_real_min>& params) const override
	{
		std::lock_guard<QMutex> _lock{m_mtxMinuit};

		// set model parameters, [ident, val, err]
		std::vector<std::tuple<std::string, std::string, std::string>> sqwparams;
		for(std::size_t paramidx=0; paramidx<params.size(); ++paramidx)
		{
			const std::string& name = std::get<0>((*m_sqwparams)[paramidx]);
			sqwparams.emplace_back(std::make_tuple(name, tl::var_to_str(params[paramidx]), ""));
		}
		m_convodlg->SetSqwParams(sqwparams);

		// start convolution simulator for new parameters and wait for it to finish
		m_convodlg->StartSim1D(true, m_seed);

		// if a stop is requested, we have no other way of getting out of here than throwing an exception...
		// if not in StartFit(), this exception will be handled at the latest by TakAppl::notify()
		if(m_convodlg->StopRequested())
			throw StopRequestedEx{"Convolution fitter stop requested."};

		return m_convodlg->GetChi2();
	}


	virtual t_real_min Up() const override
	{
		// sigma^2
		return 1.;
	}


private:
	ConvoDlg* m_convodlg = nullptr;
	const ConvoDlg::t_sqwparams* m_sqwparams = nullptr;
	unsigned int m_seed{1234};

	mutable QMutex m_mtxMinuit{};
};



/**
 * start 1d or 2d convolution fits
 */
void ConvoDlg::StartFit()
{
	if(check2dMap->isChecked())
	{
		QMessageBox::critical(this, "Error", "2D fitting is not yet implemented.");
		return;
	}


	// [ ident, type, value, error, fit? ]
	t_sqwparams sqwparams = GetSqwParams(true);
	if(sqwparams.size() == 0)
	{
		QMessageBox::critical(this, "Error", "No fit parameters defined."
			" Please set them up in the model parameters dialog (\"Parameters...\" button).");
		return;
	}


	// stop any previous fits
	Stop();
	m_atStop.store(false);


	// get fit parameters
	std::ostringstream ostrZeroErr;
	bool bAnyErrorNonZero = 0;


	std::ostringstream ostrFitParamMsg;
	ostrFitParamMsg.precision(g_iPrec);
	ostrFitParamMsg << "Using fitting parameters:\n";
	ostrFitParamMsg
		<< std::setw(15) << std::left << "Name"
		<< std::setw(15) << std::left << "Initial" 
		<< std::setw(15) << std::left << "Error" 
		<< std::setw(30) << std::left << "Limits" 
		<< "\n";

	minuit::MnUserParameters params;
	for(const auto& sqwparam : sqwparams)
	{
		// set value and error
		const std::string& varname = std::get<0>(sqwparam);
		t_real val = tl::str_to_var<t_real_min>(std::get<2>(sqwparam));
		t_real err = tl::str_to_var<t_real_min>(std::get<3>(sqwparam));
		params.Add(varname, val, err);

		// set parameter limits if given
		boost::optional<t_real_min> limLower, limUpper;
		std::string strLimits = std::get<5>(sqwparam);
		std::vector<std::string> vecLimits;
		tl::get_tokens<std::string, std::string>(strLimits, ":;|", vecLimits);
		if(vecLimits.size() == 2)
		{
			tl::trim(vecLimits[0]);
			tl::trim(vecLimits[1]);

			if(vecLimits[0]!="" && tl::str_to_lower(vecLimits[0])!="open" && tl::str_to_lower(vecLimits[0])!="none")
			{
				limLower = tl::str_to_var<t_real_min>(vecLimits[0]);
				params.SetLowerLimit(varname, *limLower);
			}
			if(vecLimits[1]!="" && tl::str_to_lower(vecLimits[1])!="open" && tl::str_to_lower(vecLimits[1])!="none")
			{
				limUpper = tl::str_to_var<t_real_min>(vecLimits[1]);
				params.SetUpperLimit(varname, *limUpper);
			}
		}

		// look for variables which have zero error (and thus cannot be fitted)
		if(tl::float_equal(err, t_real{0}))
			ostrZeroErr << varname << ", ";
		else
			bAnyErrorNonZero = 1;

		std::ostringstream ostrLim;
		ostrLim.precision(g_iPrec);
		if(limLower)
			ostrLim << *limLower;
		else
			ostrLim << "open";
		ostrLim << " : ";
		if(limUpper)
			ostrLim << *limUpper;
		else
			ostrLim << "open";

		ostrFitParamMsg
			<< std::setw(15) << std::left << varname
			<< std::setw(15) << std::left << val
			<< std::setw(15) << std::left << err
			<< std::setw(30) << std::left << ostrLim.str()
			<< "\n";
	}

	tl::log_info(ostrFitParamMsg.str());

	if(ostrZeroErr.str() != "")
	{
		std::string msg = "The error of the following parameters is zero:\n\n"
			+ ostrZeroErr.str() + "\n\nno fitting will be done for them."
			" Please set the errors to non-zero in the model parameters dialog (\"Parameters...\" button).";

		QMessageBox::warning(this, "Warning", msg.c_str());
	}

	// nothing to be done if all errors are zero
	if(!bAnyErrorNonZero)
		return;


	// minimise
	bool mini_valid = false;
	std::unique_ptr<minuit::FunctionMinimum> mini;

	try
	{
		MinuitFunc fkt{this, &sqwparams};

		minuit::MnStrategy strat(spinStrategy->value());
		std::unique_ptr<minuit::MnApplication> minimiser;

		if(comboFitter->currentIndex() == 0)
		{
			minimiser.reset(new minuit::MnSimplex(fkt, params, strat));
		}
		else if(comboFitter->currentIndex() == 1)
		{
			minimiser.reset(new minuit::MnMigrad(fkt, params, strat));
		}
		else
		{
			QMessageBox::critical(this, "Error", "Invalid minimiser.");
			return;
		}

		mini.reset(new minuit::FunctionMinimum((*minimiser)(spinMaxCalls->value(), spinTolerance->value())));
		mini_valid = mini->IsValid() && mini->HasValidParameters() && mini->UserState().IsValid();
	}
	catch(const StopRequestedEx& req)
	{
		tl::log_info(req.what());
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}


	std::ostringstream ostr_fitresults;
	ostr_fitresults.precision(g_iPrec);

	if(!mini_valid || !mini)
	{
		QMessageBox::critical(this, "Error", "Convolution fitter did not converge.");
		ostr_fitresults << "# Warning: Convolution fit did not converge.\n";
	}


	// get back minimised parameters, [ident, val, err]
	if(mini)
	{
		tl::log_debug("Final fit results:\n", *mini);
		ostr_fitresults << "# Fit chi^2: " << GetChi2() << "\n";

		std::vector<std::tuple<std::string, std::string, std::string>> newsqwparams;
		for(std::size_t paramidx=0; paramidx<sqwparams.size(); ++paramidx)
		{
			const std::string& name = std::get<0>(sqwparams[paramidx]);

			t_real_min dVal = mini->UserState().Value(name);
			t_real_min dErr = mini->UserState().Error(name);

			ostr_fitresults << "# Fitted Variable: " << name << " = " << dVal << " +- " << dErr << "\n";

			std::string val = tl::var_to_str(mini->UserState().Value(name));
			std::string err = tl::var_to_str(mini->UserState().Error(name));

			newsqwparams.emplace_back(std::make_tuple(name, val, err));
		}
		this->SetSqwParams(newsqwparams);
	}


	QString strResults = "# --------------------------------------------------------------------------------\n";
	strResults += ostr_fitresults.str().c_str();
	strResults += "# --------------------------------------------------------------------------------\n";
	strResults += textResult->toPlainText();
	textResult->setPlainText(strResults);
}
