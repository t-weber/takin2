/**
 * convolution simulation -- CLI program
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2020
 * @license GPLv2
 */

#include <boost/asio/io_service.hpp>
#include <boost/asio/signal_set.hpp>
#include <boost/scope_exit.hpp>
#include <boost/program_options.hpp>

#include <vector>
#include <thread>
#include <atomic>
#include <memory>
#include <string>

#include "monteconvo_cli.h"
#include "monteconvo_common.h"
#include "sqwfactory.h"

#include "tools/res/defs.h"
#include "tools/convofit/scan.h"
#include "TASReso.h"

#include "libs/version.h"
#include "libs/globals.h"

#include "tlibs/file/file.h"
#include "tlibs/file/prop.h"
#include "tlibs/log/log.h"
#include "tlibs/log/debug.h"
#include "tlibs/time/stopwatch.h"
#include "tlibs/helper/thread.h"
#include "tlibs/math/stat.h"
#include "tlibs/math/linalg.h"

namespace asio = boost::asio;
namespace sys = boost::system;
namespace opts = boost::program_options;

using t_real = t_real_reso;


static const std::string g_strXmlRoot("taz/");


// ----------------------------------------------------------------------------
// configuration
struct ConvoConfig
{
	t_real h_from{}, k_from{}, l_from{}, E_from{};
	t_real h_to{}, k_to{}, l_to{}, E_to{};
	t_real h_to_2{}, k_to_2{}, l_to_2{}, E_to_2{};
	t_real kfix{};
	t_real tolerance{};
	t_real S_scale{1}, S_slope{0}, S_offs{0};

	unsigned int neutron_count{500};
	unsigned int sample_step_count{1};
	unsigned int step_count{256};

	bool scan_2d{false};
	bool recycle_neutrons{true};
	bool normalise{true};
	bool flip_coords{false};
	bool has_scanfile{false};

	int algo{1};
	int fixedk{1};
	int mono_foc{1}, ana_foc{1};
	int scanaxis{4}, scanaxis2{0};

	std::string crys{}, instr{};
	std::string sqw{}, sqw_conf{};
	std::string scanfile{};
	std::string counter{}, monitor{};
	std::string temp_override{}, field_override{};
	std::string autosave{};
	std::string filter_col{}, filter_val{};
};


/**
 * loads the configuration for the convolution from a job file
 */
static ConvoConfig load_config(const tl::Prop<std::string>& xml)
{
	ConvoConfig cfg;

	// real values
	boost::optional<t_real> odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/h_from"); if(odVal) cfg.h_from = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/k_from"); if(odVal) cfg.k_from = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/l_from"); if(odVal) cfg.l_from = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/E_from"); if(odVal) cfg.E_from = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/h_to"); if(odVal) cfg.h_to = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/k_to"); if(odVal) cfg.k_to = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/l_to"); if(odVal) cfg.l_to = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/E_to"); if(odVal) cfg.E_to = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/h_to_2"); if(odVal) cfg.h_to_2 = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/k_to_2"); if(odVal) cfg.k_to_2 = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/l_to_2"); if(odVal) cfg.l_to_2 = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/E_to_2"); if(odVal) cfg.E_to_2 = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/kfix"); if(odVal) cfg.kfix = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"convofit/tolerance"); if(odVal) cfg.tolerance = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/S_scale"); if(odVal) cfg.S_scale = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/S_slope"); if(odVal) cfg.S_slope = *odVal;
	odVal = xml.QueryOpt<t_real>(g_strXmlRoot+"monteconvo/S_offs"); if(odVal) cfg.S_offs = *odVal;

	// int values
	boost::optional<int> oiVal;
	oiVal = xml.QueryOpt<unsigned int>(g_strXmlRoot+"monteconvo/neutron_count"); if(oiVal) cfg.neutron_count = *oiVal;
	oiVal = xml.QueryOpt<unsigned int>(g_strXmlRoot+"monteconvo/sample_step_count"); if(oiVal) cfg.sample_step_count = *oiVal;
	oiVal = xml.QueryOpt<unsigned int>(g_strXmlRoot+"monteconvo/step_count"); if(oiVal) cfg.step_count = *oiVal;
	//oiVal = xml.QueryOpt<unsigned int>(g_strXmlRoot+"convofit/strategy"); if(oiVal) cfg.strategy = *oiVal;
	//oiVal = xml.QueryOpt<unsigned int>(g_strXmlRoot+"convofit/max_calls"); if(oiVal) cfg.max_calls = *oiVal;

	// bool values
	boost::optional<int> obVal;
	obVal = xml.QueryOpt<int>(g_strXmlRoot+"monteconvo/scan_2d"); if(obVal) cfg.scan_2d = *obVal != 0;
	obVal = xml.QueryOpt<int>(g_strXmlRoot+"convofit/recycle_neutrons"); if(obVal) cfg.recycle_neutrons = *obVal != 0;
	obVal = xml.QueryOpt<int>(g_strXmlRoot+"convofit/normalise"); if(obVal) cfg.normalise = *obVal != 0;
	obVal = xml.QueryOpt<int>(g_strXmlRoot+"convofit/flip_coords"); if(obVal) cfg.flip_coords = *obVal != 0;
	obVal = xml.QueryOpt<int>(g_strXmlRoot+"monteconvo/has_scanfile"); if(obVal) cfg.has_scanfile = *obVal != 0;

	// indices for gui comboboxes
	boost::optional<int> oCmb;
	oCmb = xml.QueryOpt<int>(g_strXmlRoot+"monteconvo/algo"); if(oCmb) cfg.algo = *oCmb;
	oCmb = xml.QueryOpt<int>(g_strXmlRoot+"monteconvo/fixedk"); if(oCmb) cfg.fixedk = *oCmb;
	oCmb = xml.QueryOpt<int>(g_strXmlRoot+"monteconvo/mono_foc"); if(oCmb) cfg.mono_foc = *oCmb;
	oCmb = xml.QueryOpt<int>(g_strXmlRoot+"monteconvo/ana_foc"); if(oCmb) cfg.ana_foc = *oCmb;
	oCmb = xml.QueryOpt<int>(g_strXmlRoot+"convofit/scanaxis"); if(oCmb) cfg.scanaxis = *oCmb;
	oCmb = xml.QueryOpt<int>(g_strXmlRoot+"convofit/scanaxis2"); if(oCmb) cfg.scanaxis2 = *oCmb;
	//oCmb = xml.QueryOpt<int>(g_strXmlRoot+"convofit/minimiser"); if(oCmb) cfg.minimiser = *oCmb;

	// string values
	boost::optional<std::string> osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"monteconvo/crys"); if(osVal) cfg.crys = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"monteconvo/instr"); if(osVal) cfg.instr = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"monteconvo/sqw"); if(osVal) cfg.sqw = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"monteconvo/sqw_conf"); if(osVal) cfg.sqw_conf = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"monteconvo/scanfile"); if(osVal) cfg.scanfile = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"convofit/counter"); if(osVal) cfg.counter = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"convofit/monitor"); if(osVal) cfg.monitor = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"convofit/temp_override"); if(osVal) cfg.temp_override = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"convofit/field_override"); if(osVal) cfg.field_override = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"monteconvo/autosave"); if(osVal) cfg.autosave = *osVal;
	//osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"convofit/sqw_params"); if(osVal) cfg.sqw_params = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"monteconvo/filter_col"); if(osVal) cfg.filter_col = *osVal;
	osVal = xml.QueryOpt<std::string>(g_strXmlRoot+"monteconvo/filter_val"); if(osVal) cfg.filter_val = *osVal;

	return cfg;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
static std::shared_ptr<SqwBase> create_sqw_model(const std::string& strSqwIdent, const std::string& _strSqwFile)
{
	std::string strSqwFile = _strSqwFile;
	tl::trim(strSqwFile);
	strSqwFile = find_file_in_global_paths(strSqwFile);

	if(strSqwFile == "")
	{
		tl::log_warn("No S(q,w) config file given.");
		return nullptr;
	}

	std::shared_ptr<SqwBase> pSqw = construct_sqw(strSqwIdent, strSqwFile);
	if(!pSqw)
	{
		tl::log_err("Unknown S(q,w) model selected.");
		return nullptr;
	}

	if(!pSqw->IsOk())
	{
		tl::log_err("Could not create S(q,w).");
		return nullptr;
	}

	return pSqw;
}



/**
 * create 1d convolution
 */
static bool start_convo_1d(const ConvoConfig& cfg, const tl::Prop<std::string>& xml)
{
	std::shared_ptr<SqwBase> pSqw = create_sqw_model(cfg.sqw, cfg.sqw_conf);
	if(!pSqw)
		return false;
	if(!load_sqw_params(pSqw.get(), xml, g_strXmlRoot+"monteconvo/"))
		return false;

	Filter filter;
	if(cfg.filter_col != "")
		filter.colEquals = std::make_pair(cfg.filter_col, cfg.filter_val);

	Scan scan;
	if(cfg.has_scanfile)
	{
		if(!load_scan_file(cfg.scanfile, scan, cfg.flip_coords, filter))
		{
			tl::log_err("Cannot load scan(s) \"", cfg.scanfile, "\".");
			return false;
		}

		if(!scan.vecPoints.size())
		{
			tl::log_err("No points in scan(s) \"", cfg.scanfile, "\".");
			return false;
		}
	}


	std::string strAutosave = cfg.autosave;
	if(strAutosave == "")
	{
		strAutosave = "out.dat";
		tl::log_warn("Output file not set, using \"", strAutosave, "\".");
	}


	tl::Stopwatch<t_real> watch;
	watch.start();

	bool bScanAxisFound = 0;
	int iScanAxisIdx = 0;
	std::string strScanVar = "";
	std::vector<std::vector<t_real>> vecAxes;
	std::tie(bScanAxisFound, iScanAxisIdx, strScanVar, vecAxes) = get_scan_axis<t_real>(
		true, cfg.scanaxis, cfg.step_count, EPS_RLU,
		cfg.h_from, cfg.h_to, cfg.k_from, cfg.k_to, cfg.l_from, cfg.l_to, cfg.E_from, cfg.E_to);
	if(!bScanAxisFound)
	{
		tl::log_err("No scan variable found.");
		return false;
	}

	const std::vector<t_real> *pVecScanX = &vecAxes[iScanAxisIdx];
	const std::vector<t_real>& vecH = vecAxes[0];
	const std::vector<t_real>& vecK = vecAxes[1];
	const std::vector<t_real>& vecL = vecAxes[2];
	const std::vector<t_real>& vecE = vecAxes[3];


	// -------------------------------------------------------------------------
	// Load reso file
	TASReso reso;

	std::string _strResoFile = cfg.instr;
	tl::trim(_strResoFile);
	const std::string strResoFile = find_file_in_global_paths(_strResoFile);

	tl::log_debug("Loading resolution from \"", strResoFile, "\".");
	if(strResoFile == "" || !reso.LoadRes(strResoFile.c_str()))
	{
		tl::log_err("Could not load resolution file \"", strResoFile, "\".");
		return false;
	}
	// -------------------------------------------------------------------------


	if(cfg.has_scanfile)	// get crystal definition from scan file
	{
		ublas::vector<t_real> vec1 =
			tl::make_vec({scan.plane.vec1[0], scan.plane.vec1[1], scan.plane.vec1[2]});
		ublas::vector<t_real> vec2 =
			tl::make_vec({scan.plane.vec2[0], scan.plane.vec2[1], scan.plane.vec2[2]});

		reso.SetLattice(scan.sample.a, scan.sample.b, scan.sample.c,
			scan.sample.alpha, scan.sample.beta, scan.sample.gamma,
			vec1, vec2);
	}
	else	// use crystal config file
	{
		// -------------------------------------------------------------------------
		// Load lattice
		std::string _strLatticeFile = cfg.crys;
		tl::trim(_strLatticeFile);
		const std::string strLatticeFile = find_file_in_global_paths(_strLatticeFile);

		tl::log_debug("Loading crystal from \"", strLatticeFile, "\".");
		if(strLatticeFile == "" || !reso.LoadLattice(strLatticeFile.c_str()))
		{
			tl::log_err("Could not load crystal file \"", strLatticeFile, "\".");
			return false;
		}
		// -------------------------------------------------------------------------
	}

	reso.SetAlgo(ResoAlgo(cfg.algo+1));
	reso.SetKiFix(cfg.fixedk==0);
	reso.SetKFix(cfg.kfix);
	reso.SetOptimalFocus(get_reso_focus(cfg.mono_foc, cfg.ana_foc));


	std::ostringstream ostrOut;
	ostrOut.precision(g_iPrec);
	ostrOut << "#\n";
	ostrOut << "# Takin/Monteconvo version " << TAKIN_VER << "\n";
	ostrOut << "# MC neutrons: " << cfg.neutron_count << "\n";
	ostrOut << "# MC sample steps: " << cfg.sample_step_count << "\n";
	ostrOut << "# Scale: " << cfg.S_scale << "\n";
	ostrOut << "# Slope: " << cfg.S_slope << "\n";
	ostrOut << "# Offset: " << cfg.S_offs << "\n";
	if(cfg.scanfile != "")
		ostrOut << "# Scan file: " << cfg.scanfile << "\n";
	ostrOut << "#\n";

	ostrOut << std::left << std::setw(g_iPrec*2) << "# h" << " "
		<< std::left << std::setw(g_iPrec*2) << "k" << " "
		<< std::left << std::setw(g_iPrec*2) << "l" << " "
		<< std::left << std::setw(g_iPrec*2) << "E" << " "
		<< std::left << std::setw(g_iPrec*2) << "S(Q,E)" << "\n";


	std::vector<t_real_reso> vecQ, vecS, vecScaledS;

	vecQ.reserve(cfg.step_count);
	vecS.reserve(cfg.step_count);
	vecScaledS.reserve(cfg.step_count);

	unsigned int iNumThreads = get_max_threads();
	tl::log_debug("Calculating using ", iNumThreads, " threads.");

	void (*pThStartFunc)() = []{ tl::init_rand(); };
	tl::ThreadPool<std::pair<bool, t_real>()> tp(iNumThreads, pThStartFunc);
	auto& lstFuts = tp.GetResults();

	for(unsigned int iStep=0; iStep<cfg.step_count; ++iStep)
	{
		t_real dCurH = vecH[iStep];
		t_real dCurK = vecK[iStep];
		t_real dCurL = vecL[iStep];
		t_real dCurE = vecE[iStep];

		tp.AddTask([&reso, dCurH, dCurK, dCurL, dCurE, pSqw, &cfg]()
			-> std::pair<bool, t_real>
		{
			t_real dS = 0.;
			t_real dhklE_mean[4] = {0., 0., 0., 0.};

			if(cfg.neutron_count == 0)
			{	// if no neutrons are given, just plot the unconvoluted S(q,w)
				// TODO: add an option to let the user choose if S(Q,E) is
				// really the dynamical structure factor, or its absolute square
				dS += (*pSqw)(dCurH, dCurK, dCurL, dCurE);
			}
			else
			{	// convolution
				TASReso localreso = reso;
				localreso.SetRandomSamplePos(cfg.sample_step_count);
				std::vector<ublas::vector<t_real>> vecNeutrons;

				try
				{
					if(!localreso.SetHKLE(dCurH, dCurK, dCurL, dCurE))
					{
						std::ostringstream ostrErr;
						ostrErr << "Invalid crystal position: (" <<
							dCurH << " " << dCurK << " " << dCurL << ") rlu, "
							<< dCurE << " meV.";
						throw tl::Err(ostrErr.str().c_str());
					}
				}
				catch(const std::exception& ex)
				{
					tl::log_err(ex.what());
					return std::pair<bool, t_real>(false, 0.);
				}

				Ellipsoid4d<t_real> elli =
					localreso.GenerateMC_deferred(cfg.neutron_count, vecNeutrons);

				for(const ublas::vector<t_real>& vecHKLE : vecNeutrons)
				{
					// TODO: add an option to let the user choose if S(Q,E) is
					// really the dynamical structure factor, or its absolute square
					dS += (*pSqw)(vecHKLE[0], vecHKLE[1], vecHKLE[2], vecHKLE[3]);

					for(int i=0; i<4; ++i)
						dhklE_mean[i] += vecHKLE[i];
				}

				dS /= t_real(cfg.neutron_count*cfg.sample_step_count);
				for(int i=0; i<4; ++i)
					dhklE_mean[i] /= t_real(cfg.neutron_count*cfg.sample_step_count);

				if(localreso.GetResoParams().flags & CALC_R0)
					dS *= localreso.GetResoResults().dR0;
				if(localreso.GetResoParams().flags & CALC_RESVOL)
					dS /= localreso.GetResoResults().dResVol * tl::get_pi<t_real>() * t_real(3.);
			}
			return std::pair<bool, t_real>(true, dS);
		});
	}

	tp.Start();
	auto iterTask = tp.GetTasks().begin();
	unsigned int iStep = 0;
	for(auto &fut : lstFuts)
	{
		// deferred (in main thread), eval this task manually
		if(iNumThreads == 0)
		{
			(*iterTask)();
			++iterTask;
		}

		std::pair<bool, t_real> pairS = fut.get();
		if(!pairS.first) break;
		t_real dS = pairS.second;
		if(tl::is_nan_or_inf(dS))
		{
			dS = t_real(0);
			tl::log_warn("S(q,w) is invalid.");
		}

		ostrOut << std::left << std::setw(g_iPrec*2) << vecH[iStep] << " "
			<< std::left << std::setw(g_iPrec*2) << vecK[iStep] << " "
			<< std::left << std::setw(g_iPrec*2) << vecL[iStep] << " "
			<< std::left << std::setw(g_iPrec*2) << vecE[iStep] << " "
			<< std::left << std::setw(g_iPrec*2) << dS << "\n";

		const t_real dXVal = (*pVecScanX)[iStep];
		t_real dYVal = cfg.S_scale*(dS + cfg.S_slope*dXVal) + cfg.S_offs;
		if(dYVal < 0.)
			dYVal = 0.;

		vecQ.push_back(dXVal);
		vecS.push_back(dS);
		vecScaledS.push_back(dYVal);

		static const std::vector<t_real> vecNull;
		bool bIsLastStep = (iStep == lstFuts.size()-1);

		if(bIsLastStep)
		{
			if(bIsLastStep)
				ostrOut << "# ------------------------- EOF -------------------------\n";

			// output
			std::ofstream ofstrAutosave(strAutosave);
			ofstrAutosave << ostrOut.str() << std::endl;
		}


		std::string strStopTime = watch.GetEstStopTimeStr(t_real(iStep+1)/t_real(cfg.step_count));
		std::cout << "\rStep " << iStep+1 << "/" << cfg.step_count << ". Estimated stop time: " << strStopTime << "...          ";
		if(iStep+1 == cfg.step_count)
			std::cout << "\n";
		std::cout.flush();

		++iStep;
	}
	tl::log_info("Convolution simulation finished.");


	// approximate chi^2
	if(cfg.has_scanfile && pSqw)
	{
		const std::size_t iNumScanPts = scan.vecPoints.size();
		std::vector<t_real> vecSFuncY;
		vecSFuncY.reserve(iNumScanPts);

		for(std::size_t iScanPt=0; iScanPt<iNumScanPts; ++iScanPt)
		{
			const ScanPoint& pt = scan.vecPoints[iScanPt];
			t_real E = pt.E / tl::one_meV;
			ublas::vector<t_real> vecScanHKLE = tl::make_vec({ pt.h, pt.k, pt.l, E });


			// find point on S(q,w) curve closest to scan point
			std::size_t iMinIdx = 0;
			t_real dMinDist = std::numeric_limits<t_real>::max();
			for(std::size_t iStep=0; iStep<cfg.step_count; ++iStep)
			{
				ublas::vector<t_real> vecCurveHKLE =
					tl::make_vec({ vecH[iStep], vecK[iStep], vecL[iStep], vecE[iStep] });

				t_real dDist = ublas::norm_2(vecCurveHKLE - vecScanHKLE);
				if(dDist < dMinDist)
				{
					dMinDist = dDist;
					iMinIdx = iStep;
				}
			}

			// add the scaled S value from the closest point
			vecSFuncY.push_back(vecScaledS[iMinIdx]);
		}

		t_real tChi2 = tl::chi2_direct<t_real>(iNumScanPts,
			vecSFuncY.data(), scan.vecCts.data(), scan.vecCtsErr.data());
		tl::log_info("chi^2 = ", tChi2);

		// output
		std::ofstream ofstrAutosave(strAutosave, std::ios_base::app);
		ofstrAutosave << "# chi^2: " << tChi2 << std::endl;
	}


	// output elapsed time
	watch.stop();

	// output
	std::ofstream ofstrAutosave(strAutosave, std::ios_base::app);
	ofstrAutosave << "# Simulation start time: " << watch.GetStartTimeStr() << "\n";
	ofstrAutosave << "# Simulation stop time: " << watch.GetStopTimeStr() << std::endl;

	return true;
}



/**
 * create 2d convolution
 */
static bool start_convo_2d(const ConvoConfig& cfg, const tl::Prop<std::string>& xml)
{
	std::shared_ptr<SqwBase> pSqw = create_sqw_model(cfg.sqw, cfg.sqw_conf);
	if(!pSqw)
		return false;
	if(!load_sqw_params(pSqw.get(), xml, g_strXmlRoot+"monteconvo/"))
		return false;


	std::string strAutosave = cfg.autosave;
	if(strAutosave == "")
	{
		strAutosave = "out.dat";
		tl::log_warn("Output file not set, using \"", strAutosave, "\".");
	}


	tl::Stopwatch<t_real> watch;
	watch.start();

	const t_real dStartHKL[] = { cfg.h_from, cfg.k_from, cfg.l_from, cfg.E_from };
	const t_real dDeltaHKL1[] =
	{
		(cfg.h_to - cfg.h_from) / t_real(cfg.step_count),
		(cfg.k_to - cfg.k_from) / t_real(cfg.step_count),
		(cfg.l_to - cfg.l_from) / t_real(cfg.step_count),
		(cfg.E_to - cfg.E_from) / t_real(cfg.step_count)
	};
	const t_real dDeltaHKL2[] =
	{
		(cfg.h_to_2 - cfg.h_from) / t_real(cfg.step_count),
		(cfg.k_to_2 - cfg.k_from) / t_real(cfg.step_count),
		(cfg.l_to_2 - cfg.l_from) / t_real(cfg.step_count),
		(cfg.E_to_2 - cfg.E_from) / t_real(cfg.step_count)
	};


	// -------------------------------------------------------------------------
	// find axis labels and ranges
	std::string strScanVar1 = "";
	t_real dStart1{}, dStop1{};
	if(cfg.scanaxis==1 || (cfg.scanaxis==0 && !tl::float_equal(cfg.h_from, cfg.h_to, EPS_RLU)))
	{
		strScanVar1 = "h (rlu)";
		dStart1 = cfg.h_from;
		dStop1 = cfg.h_to;
	}
	else if(cfg.scanaxis==2 || (cfg.scanaxis==0 && !tl::float_equal(cfg.k_from, cfg.k_to, EPS_RLU)))
	{
		strScanVar1 = "k (rlu)";
		dStart1 = cfg.k_from;
		dStop1 = cfg.k_to;
	}
	else if(cfg.scanaxis==3 || (cfg.scanaxis==0 && !tl::float_equal(cfg.l_from, cfg.l_to, EPS_RLU)))
	{
		strScanVar1 = "l (rlu)";
		dStart1 = cfg.l_from;
		dStop1 = cfg.l_to;
	}
	else if(cfg.scanaxis==4 || (cfg.scanaxis==0 && !tl::float_equal(cfg.E_from, cfg.E_to, EPS_RLU)))
	{
		strScanVar1 = "E (meV)";
		dStart1 = cfg.E_from;
		dStop1 = cfg.E_to;
	}

	std::string strScanVar2 = "";
	t_real dStart2{}, dStop2{};
	if(cfg.scanaxis2==1 || (cfg.scanaxis2==0 && !tl::float_equal(cfg.h_from, cfg.h_to_2, EPS_RLU)))
	{
		strScanVar2 = "h (rlu)";
		dStart2 = cfg.h_from;
		dStop2 = cfg.h_to_2;
	}
	else if(cfg.scanaxis2==2 || (cfg.scanaxis2==0 && !tl::float_equal(cfg.k_from, cfg.k_to_2, EPS_RLU)))
	{
		strScanVar2 = "k (rlu)";
		dStart2 = cfg.k_from;
		dStop2 = cfg.k_to_2;
	}
	else if(cfg.scanaxis2==3 || (cfg.scanaxis2==0 && !tl::float_equal(cfg.l_from, cfg.l_to_2, EPS_RLU)))
	{
		strScanVar2 = "l (rlu)";
		dStart2 = cfg.l_from;
		dStop2 = cfg.l_to_2;
	}
	else if(cfg.scanaxis2==4 || (cfg.scanaxis2==0 && !tl::float_equal(cfg.E_from, cfg.E_to_2, EPS_RLU)))
	{
		strScanVar2 = "E (meV)";
		dStart2 = cfg.E_from;
		dStop2 = cfg.E_to_2;
	}

	// -------------------------------------------------------------------------



	// -------------------------------------------------------------------------
	// Load reso file
	TASReso reso;
	std::string _strResoFile = cfg.instr;
	tl::trim(_strResoFile);
	const std::string strResoFile = find_file_in_global_paths(_strResoFile);

	tl::log_debug("Loading resolution from \"", strResoFile, "\".");
	if(strResoFile == "" || !reso.LoadRes(strResoFile.c_str()))
	{
		tl::log_err("Could not load resolution file \"", strResoFile, "\".");
		return false;
	}
	// -------------------------------------------------------------------------


	// -------------------------------------------------------------------------
	// Load lattice
	std::string _strLatticeFile = cfg.crys;
	tl::trim(_strLatticeFile);
	const std::string strLatticeFile = find_file_in_global_paths(_strLatticeFile);

	tl::log_debug("Loading crystal from \"", strLatticeFile, "\".");
	if(strLatticeFile == "" || !reso.LoadLattice(strLatticeFile.c_str()))
	{
		tl::log_err("Could not load crystal file \"", strLatticeFile, "\".");
		return false;
	}
	// -------------------------------------------------------------------------


	reso.SetAlgo(ResoAlgo(cfg.algo+1));
	reso.SetKiFix(cfg.fixedk==0);
	reso.SetKFix(cfg.kfix);
	reso.SetOptimalFocus(get_reso_focus(cfg.mono_foc, cfg.ana_foc));


	std::ostringstream ostrOut;
	ostrOut.precision(g_iPrec);
	ostrOut << "#\n";
	ostrOut << "# Takin/Monteconvo version " << TAKIN_VER << "\n";
	ostrOut << "# MC neutrons: " << cfg.neutron_count << "\n";
	ostrOut << "# MC sample steps: " << cfg.sample_step_count << "\n";
	ostrOut << "#\n";
	ostrOut << std::left << std::setw(g_iPrec*2) << "# h" << " "
		<< std::left << std::setw(g_iPrec*2) << "k" << " "
		<< std::left << std::setw(g_iPrec*2) << "l" << " "
		<< std::left << std::setw(g_iPrec*2) << "E" << " "
		<< std::left << std::setw(g_iPrec*2) << "S(Q,E)" << "\n";


	std::vector<t_real> vecH; vecH.reserve(cfg.step_count*cfg.step_count);
	std::vector<t_real> vecK; vecK.reserve(cfg.step_count*cfg.step_count);
	std::vector<t_real> vecL; vecL.reserve(cfg.step_count*cfg.step_count);
	std::vector<t_real> vecE; vecE.reserve(cfg.step_count*cfg.step_count);

	for(unsigned int iStepY=0; iStepY<cfg.step_count; ++iStepY)
	{
		for(unsigned int iStepX=0; iStepX<cfg.step_count; ++iStepX)
		{
			vecH.push_back(dStartHKL[0] + dDeltaHKL2[0]*t_real(iStepY) + dDeltaHKL1[0]*t_real(iStepX));
			vecK.push_back(dStartHKL[1] + dDeltaHKL2[1]*t_real(iStepY) + dDeltaHKL1[1]*t_real(iStepX));
			vecL.push_back(dStartHKL[2] + dDeltaHKL2[2]*t_real(iStepY) + dDeltaHKL1[2]*t_real(iStepX));
			vecE.push_back(dStartHKL[3] + dDeltaHKL2[3]*t_real(iStepY) + dDeltaHKL1[3]*t_real(iStepX));
		}
	}

	unsigned int iNumThreads = get_max_threads();
	tl::log_debug("Calculating using ", iNumThreads, " threads.");

	void (*pThStartFunc)() = []{ tl::init_rand(); };
	tl::ThreadPool<std::pair<bool, t_real>()> tp(iNumThreads, pThStartFunc);
	auto& lstFuts = tp.GetResults();

	for(unsigned int iStep=0; iStep<cfg.step_count*cfg.step_count; ++iStep)
	{
		t_real dCurH = vecH[iStep];
		t_real dCurK = vecK[iStep];
		t_real dCurL = vecL[iStep];
		t_real dCurE = vecE[iStep];

		tp.AddTask(
		[&reso, dCurH, dCurK, dCurL, dCurE, pSqw, &cfg]()
			-> std::pair<bool, t_real>
		{
			t_real dS = 0.;
			t_real dhklE_mean[4] = {0., 0., 0., 0.};

			if(cfg.neutron_count == 0)
			{	// if no neutrons are given, just plot the unconvoluted S(q,w)
				// TODO: add an option to let the user choose if S(Q,E) is
				// really the dynamical structure factor, or its absolute square
				dS += (*pSqw)(dCurH, dCurK, dCurL, dCurE);
			}
			else
			{	// convolution
				TASReso localreso = reso;
				localreso.SetRandomSamplePos(cfg.sample_step_count);
				std::vector<ublas::vector<t_real>> vecNeutrons;

				try
				{
					if(!localreso.SetHKLE(dCurH, dCurK, dCurL, dCurE))
					{
						std::ostringstream ostrErr;
						ostrErr << "Invalid crystal position: (" <<
							dCurH << " " << dCurK << " " << dCurL << ") rlu, "
							<< dCurE << " meV.";
						throw tl::Err(ostrErr.str().c_str());
					}
				}
				catch(const std::exception& ex)
				{
					//QMessageBox::critical(this, "Error", ex.what());
					tl::log_err(ex.what());
					return std::pair<bool, t_real>(false, 0.);
				}

				Ellipsoid4d<t_real> elli =
					localreso.GenerateMC_deferred(cfg.neutron_count, vecNeutrons);

				for(const ublas::vector<t_real>& vecHKLE : vecNeutrons)
				{
					// TODO: add an option to let the user choose if S(Q,E) is
					// really the dynamical structure factor, or its absolute square
					dS += (*pSqw)(vecHKLE[0], vecHKLE[1], vecHKLE[2], vecHKLE[3]);

					for(int i=0; i<4; ++i)
						dhklE_mean[i] += vecHKLE[i];
				}

				dS /= t_real(cfg.neutron_count*cfg.sample_step_count);
				for(int i=0; i<4; ++i)
					dhklE_mean[i] /= t_real(cfg.neutron_count*cfg.sample_step_count);

				if(localreso.GetResoParams().flags & CALC_R0)
					dS *= localreso.GetResoResults().dR0;
				if(localreso.GetResoParams().flags & CALC_RESVOL)
					dS /= localreso.GetResoResults().dResVol * tl::get_pi<t_real>() * t_real(3.);
			}
			return std::pair<bool, t_real>(true, dS);
		});
	}

	tp.Start();
	auto iterTask = tp.GetTasks().begin();
	unsigned int iStep = 0;
	for(auto &fut : lstFuts)
	{
		// deferred (in main thread), eval this task manually
		if(iNumThreads == 0)
		{
			(*iterTask)();
			++iterTask;
		}

		std::pair<bool, t_real> pairS = fut.get();
		if(!pairS.first) break;
		t_real dS = pairS.second;
		if(tl::is_nan_or_inf(dS))
		{
			dS = t_real(0);
			tl::log_warn("S(q,w) is invalid.");
		}

		ostrOut << std::left << std::setw(g_iPrec*2) << vecH[iStep] << " "
			<< std::left << std::setw(g_iPrec*2) << vecK[iStep] << " "
			<< std::left << std::setw(g_iPrec*2) << vecL[iStep] << " "
			<< std::left << std::setw(g_iPrec*2) << vecE[iStep] << " "
			<< std::left << std::setw(g_iPrec*2) << dS << "\n";


		bool bIsLastStep = (iStep == lstFuts.size()-1);

		if(bIsLastStep)
		{
			if(bIsLastStep)
				ostrOut << "# ------------------------- EOF -------------------------\n";

			// output
			std::ofstream ofstrAutosave(strAutosave);
			ofstrAutosave << ostrOut.str() << std::endl;
		}

		std::string strStopTime = watch.GetEstStopTimeStr(t_real(iStep+1)/t_real(cfg.step_count*cfg.step_count));
		std::cout << "\rStep " << iStep+1 << "/" << cfg.step_count << ". Estimated stop time: " << strStopTime << "...          ";
		if(iStep+1 == cfg.step_count*cfg.step_count)
			std::cout << "\n";
		std::cout.flush();

		++iStep;
	}
	tl::log_info("Convolution simulation finished.");

	// output elapsed time
	watch.stop();

	// output
	std::ofstream ofstrAutosave(strAutosave, std::ios_base::app);
	ofstrAutosave << "# Simulation start time: " << watch.GetStartTimeStr() << "\n";
	ofstrAutosave << "# Simulation stop time: " << watch.GetStopTimeStr() << std::endl;

	return true;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// main program

int monteconvo_main(int argc, char** argv)
{
	try
	{
#ifdef MONTECONVO_STANDALONE	// only show copyright banner if not already displayed from Takin main program
		tl::log_info("--------------------------------------------------------------------------------");
		tl::log_info("This is the Takin command-line convolution simulator (monteconvo), version " TAKIN_VER ".");
		tl::log_info("Written by Tobias Weber <tweber@ill.fr>, 2014 - 2021.");
		tl::log_info(TAKIN_LICENSE("Takin/Monteconvo"));
		tl::log_debug("Resolution calculation uses ", sizeof(t_real_reso)*8, " bit ", tl::get_typename<t_real_reso>(), "s.");
		tl::log_info("--------------------------------------------------------------------------------");
#endif

		load_sqw_plugins();


		// --------------------------------------------------------------------
		// get job files and program options
		std::vector<std::string> vecJobs;

		// normal args
        opts::options_description args("monteconvo options (overriding config file settings)");
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("job-file",
			opts::value<decltype(vecJobs)>(&vecJobs),
			"convolution config file")));
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("max-threads",
			opts::value<decltype(g_iMaxThreads)>(&g_iMaxThreads),
			"maximum number of threads")));
		// dummy arg if launched from takin executable
		bool bStartedFromTakin = 0;
#ifndef MONTECONVO_STANDALONE
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("convosim",
			opts::bool_switch(&bStartedFromTakin),
			"launch convofit from takin")));
#endif


		// positional args
		opts::positional_options_description args_pos;
		args_pos.add("job-file", -1);

		opts::basic_command_line_parser<char> clparser(argc, argv);
		clparser.options(args);
		clparser.positional(args_pos);
		opts::basic_parsed_options<char> parsedopts = clparser.run();

		opts::variables_map opts_map;
		opts::store(parsedopts, opts_map);
		opts::notify(opts_map);


		int args_to_ignore = 1;	// started with "monteconvo-cli"
		if(bStartedFromTakin)
			++args_to_ignore;	// started with "takin --monteconvo-cli"

		if(argc <= args_to_ignore)
		{
			std::ostringstream ostrHelp;
			ostrHelp << "Usage: ";
			for(int argidx=0; argidx<args_to_ignore; ++argidx)
				ostrHelp << argv[argidx] << " ";
			ostrHelp << "[options] <config file>\n";
			ostrHelp << args;
			tl::log_info(ostrHelp.str());
			return -1;
		}

		if(vecJobs.size() == 0)
		{
			tl::log_err("No config files given.");
			return -1;
		}
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// load convolution job file
		const std::string& strJobFile = vecJobs[0];
		if(!tl::file_exists(strJobFile.c_str()))
		{
			tl::log_err("Convolution config file \"", strJobFile, "\" does not exist.");
			return -1;
		}

		// add the location of the convo file as a possible global path
		std::string strGlobDir = tl::get_dir(strJobFile);
		clear_global_paths();
		if(strGlobDir != "")
			add_global_path(strGlobDir);


		tl::Prop<std::string> xml;
		if(!xml.Load(strJobFile, tl::PropType::XML))
		{
			tl::log_err("Convolution config file \"", strJobFile, "\" could not be loaded.");
			return -1;
		}


		ConvoConfig cfg = load_config(xml);
		// --------------------------------------------------------------------



		tl::Stopwatch<t_real> watch;
		watch.start();

		bool ok = 0;
		if(cfg.scan_2d)
		{
			tl::log_info("Performing a 2d convolution simulation.");
			ok = start_convo_2d(cfg, xml);
		}
		else
		{
			tl::log_info("Performing a 1d convolution simulation.");
			ok = start_convo_1d(cfg, xml);
		}

		if(!ok)
			tl::log_err("Simulation failed!");


		watch.stop();
		tl::log_info("================================================================================");
		tl::log_info("Start time:     ", watch.GetStartTimeStr());
		tl::log_info("Stop time:      ", watch.GetStopTimeStr());
		tl::log_info("Execution time: ", tl::get_duration_str_secs<t_real>(watch.GetDur()));
		tl::log_info("================================================================================");

		return ok ? 0 : -1;
    }
	catch(const std::exception& ex)
	{
		tl::log_crit(ex.what());
	}

    return 0;
}
// ----------------------------------------------------------------------------
