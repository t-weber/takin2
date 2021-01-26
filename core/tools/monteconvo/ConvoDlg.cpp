/**
 * monte carlo convolution tool
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015-2020
 * @license GPLv2
 */

#include "ConvoDlg.h"
#include "tlibs/string/string.h"
#include "tlibs/math/math.h"
#include "tlibs/math/rand.h"
#include "tlibs/file/file.h"

#include "libs/globals.h"
#include "libs/globals_qt.h"
#include "libs/qt/recent.h"

#include <boost/predef.h>

#include <iostream>
#include <fstream>
#include <tuple>

#include <QFileDialog>
#include <QMessageBox>
#include <QMenu>


using t_real = t_real_reso;
const std::string ConvoDlg::s_strTitle = "Resolution Convolution";


ConvoDlg::ConvoDlg(QWidget* pParent, QSettings* pSett)
	: QDialog(pParent, Qt::WindowTitleHint|Qt::WindowCloseButtonHint|Qt::WindowMinMaxButtonsHint),
		m_pSett(pSett)
{
	setupUi(this);
	setWindowTitle(s_strTitle.c_str());

	// -------------------------------------------------------------------------
	// widgets
	m_vecSpinBoxes = { spinStartH, spinStartK, spinStartL, spinStartE,
		spinStopH, spinStopK, spinStopL, spinStopE,
		spinStopH2, spinStopK2, spinStopL2, spinStopE2,
		spinKfix,
		spinTolerance
	};

	m_vecSpinNames = {
		"monteconvo/h_from", "monteconvo/k_from", "monteconvo/l_from", "monteconvo/E_from",
		"monteconvo/h_to", "monteconvo/k_to", "monteconvo/l_to", "monteconvo/E_to",
		"monteconvo/h_to_2", "monteconvo/k_to_2", "monteconvo/l_to_2", "monteconvo/E_to_2",
		"monteconvo/kfix",
		"convofit/tolerance"
	};

	m_vecIntSpinBoxes = { spinNeutrons, spinSampleSteps,
		spinStepCnt,
		spinStrategy, spinMaxCalls
	};
	m_vecIntSpinNames = {
		"monteconvo/neutron_count", "monteconvo/sample_step_count",
		"monteconvo/step_count",
		"convofit/strategy", "convofit/max_calls"
	};

	m_vecEditBoxes = {
		editFilterCol, editFilterVal,
		editCrys, editRes, editSqw, editScan,
		editScale, editSlope, editOffs,
		editCounter, editMonitor,
		editTemp, editField,
		editAutosave,
	};
	m_vecEditNames = {
		"monteconvo/filter_col", "monteconvo/filter_val",
		"monteconvo/crys", "monteconvo/instr",
		"monteconvo/sqw_conf", "monteconvo/scanfile",
		"monteconvo/S_scale", "monteconvo/S_slope", "monteconvo/S_offs",
		"convofit/counter", "convofit/monitor",
		"convofit/temp_override", "convofit/field_override",
		"monteconvo/autosave",
	};

	m_vecTextBoxes = { editSqwParams };
	m_vecTextNames = { "convofit/sqw_params" };

	m_vecComboBoxes = { comboAlgo, comboFixedK, comboFocMono, comboFocAna,
		comboFitter, comboAxis, comboAxis2,
	};
	m_vecComboNames = { "monteconvo/algo", "monteconvo/fixedk", "monteconvo/mono_foc",
		"monteconvo/ana_foc", "convofit/minimiser", "convofit/scanaxis", "convofit/scanaxis2",
	};

	m_vecCheckBoxes = { checkScan, check2dMap,
		checkRnd, checkNorm, checkFlip
	};
	m_vecCheckNames = { "monteconvo/has_scanfile", "monteconvo/scan_2d",
		"convofit/recycle_neutrons", "convofit/normalise", "convofit/flip_coords"
	};
	// -------------------------------------------------------------------------

	if(m_pSett)
	{
		QFont font;
		if(m_pSett->contains("main/font_gen") && font.fromString(m_pSett->value("main/font_gen", "").toString()))
			setFont(font);
	}

	btnStart->setIcon(load_icon("res/icons/media-playback-start.svg"));
	btnStartFit->setIcon(load_icon("res/icons/media-playback-start.svg"));
	btnStop->setIcon(load_icon("res/icons/media-playback-stop.svg"));

	/*
	 * curve 0,1	->	convolution
	 * curve 2	->	scan points
	 * curve 3-15	->	dispersion branches
	 */
	m_plotwrap.reset(new QwtPlotWrapper(plot, CONVO_MAX_CURVES, true));
	m_plotwrap->GetPlot()->setAxisTitle(QwtPlot::xBottom, "");
	m_plotwrap->GetPlot()->setAxisTitle(QwtPlot::yLeft, "S(Q,E) (a.u.)");

	m_plotwrap2d.reset(new QwtPlotWrapper(plot2d, 1, 0, 0, 1));
	m_plotwrap2d->GetPlot()->setAxisTitle(QwtPlot::yRight, "S(Q,E) (a.u.)");


	// --------------------------------------------------------------------
	// convolution lines
	QPen penCurve;
	penCurve.setColor(QColor(0,0,0x99));
	penCurve.setWidth(2);
	m_plotwrap->GetCurve(0)->setPen(penCurve);
	m_plotwrap->GetCurve(0)->setStyle(QwtPlotCurve::CurveStyle::Lines);
	m_plotwrap->GetCurve(0)->setTitle("S(Q,E)");

	// convolution points
	QPen penPoints;
	penPoints.setColor(QColor(0xff,0,0));
	penPoints.setWidth(4);
	m_plotwrap->GetCurve(1)->setPen(penPoints);
	m_plotwrap->GetCurve(1)->setStyle(QwtPlotCurve::CurveStyle::Dots);
	m_plotwrap->GetCurve(1)->setTitle("S(Q,E)");

	// scan data points
	QPen penScanPoints;
	penScanPoints.setColor(QColor(0x00,0x90,0x00));
	penScanPoints.setWidth(6);
	m_plotwrap->GetCurve(2)->setPen(penScanPoints);
	m_plotwrap->GetCurve(2)->setStyle(QwtPlotCurve::CurveStyle::Dots);
	m_plotwrap->GetCurve(2)->setTitle("S(Q,E)");
	m_plotwrap->GetCurve(2)->SetShowErrors(true);

	// dispersion branches
	for(int iCurve=CONVO_DISP_CURVE_START; iCurve<CONVO_MAX_CURVES; ++iCurve)
	{
		m_plotwrap->GetCurve(iCurve)->setPen(penCurve);
		m_plotwrap->GetCurve(iCurve)->setStyle(QwtPlotCurve::CurveStyle::Lines);
		m_plotwrap->GetCurve(iCurve)->setTitle("E(Q)");
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	QMenu *pMenuActions = new QMenu("Actions", this);

	QAction *pHK = new QAction("h <-> k", pMenuActions);
	QAction *pHL = new QAction("h <-> l", pMenuActions);
	QAction *pKL = new QAction("k <-> l", pMenuActions);
	pMenuActions->addAction(pHK);
	pMenuActions->addAction(pHL);
	pMenuActions->addAction(pKL);

	btnActions->setMenu(pMenuActions);
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// fill sqw combo box
	load_sqw_plugins();
	auto vecSqwNames = get_sqw_names();
	for(const auto& tupSqw : vecSqwNames)
	{
		QString strIdent = std::get<0>(tupSqw).c_str();
		QString strName = std::get<1>(tupSqw).c_str();

		comboSqw->addItem(strName, strIdent);
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// menu bar
	m_pMenuBar = new QMenuBar(this);
	if(m_pSett)
		m_pMenuBar->setNativeMenuBar(m_pSett->value("main/native_dialogs", 1).toBool());


	// file menu
	QMenu *pMenuFile = new QMenu("File", this);

	QAction *pNew = new QAction("New", this);
	pNew->setIcon(load_icon("res/icons/document-new.svg"));
	pMenuFile->addAction(pNew);

	pMenuFile->addSeparator();

	QAction *pLoad = new QAction("Load...", this);
	pLoad->setIcon(load_icon("res/icons/document-open.svg"));
	pMenuFile->addAction(pLoad);

	m_pMenuRecent = new QMenu("Recently Loaded", this);
	RecentFiles recent(m_pSett, "monteconvo/recent");
	recent.FillMenu(m_pMenuRecent, [this](const std::string& str){ Load(str.c_str()); });
	pMenuFile->addMenu(m_pMenuRecent);

	QAction *pSave = new QAction("Save", this);
	pSave->setIcon(load_icon("res/icons/document-save.svg"));
	pMenuFile->addAction(pSave);

	QAction *pSaveAs = new QAction("Save As...", this);
	pSaveAs->setIcon(load_icon("res/icons/document-save-as.svg"));
	pMenuFile->addAction(pSaveAs);

	pMenuFile->addSeparator();

	QAction *pConvofit = new QAction("Export to Convofit...", this);
	pConvofit->setIcon(load_icon("res/icons/drive-harddisk.svg"));
	pMenuFile->addAction(pConvofit);

	pMenuFile->addSeparator();

	QAction *pExit = new QAction("Quit Convo", this);
	pExit->setIcon(load_icon("res/icons/system-log-out.svg"));
	pMenuFile->addAction(pExit);


	// view menu
	QMenu *pMenuView = new QMenu("View", this);

	QAction *pActionLogY = new QAction("Toggle Logarithmic Scale", this);
	pMenuView->addAction(pActionLogY);


	// actions menu
	QMenu *pMenuConvoActions = new QMenu("Actions", this);

	QAction *pActionStart = new QAction("Start Convolution Simulation", this);
	pActionStart->setIcon(load_icon("res/icons/media-playback-start.svg"));
	pMenuConvoActions->addAction(pActionStart);

	QAction *pActionStartFit = new QAction("Start Convolution Fit", this);
	pActionStartFit->setIcon(load_icon("res/icons/media-playback-start.svg"));
	pMenuConvoActions->addAction(pActionStartFit);

	QAction *pActionDisp = new QAction("Plot Dispersion", this);
	pMenuConvoActions->addAction(pActionDisp);


	// results menu
	QMenu *pMenuPlots = new QMenu("Results", this);

	m_pLiveResults = new QAction("Live Results", this);
	m_pLiveResults->setCheckable(1);
	m_pLiveResults->setChecked(0);
	pMenuPlots->addAction(m_pLiveResults);

	m_pLivePlots = new QAction("Live Plots", this);
	m_pLivePlots->setCheckable(1);
	m_pLivePlots->setChecked(1);
	pMenuPlots->addAction(m_pLivePlots);

	pMenuPlots->addSeparator();

	QAction *pExportPlot = new QAction("Export Plot Data...", this);
	pMenuPlots->addAction(pExportPlot);

	QAction *pExportPlotGpl = new QAction("Export Plot to Gnuplot...", this);
	pMenuPlots->addAction(pExportPlotGpl);

	pMenuPlots->addSeparator();

	QAction *pExportPlot2d = new QAction("Export 2D Plot Data...", this);
	pMenuPlots->addAction(pExportPlot2d);

	QAction *pExportPlot2dGpl = new QAction("Export 2D Plot to Gnuplot...", this);
	pMenuPlots->addAction(pExportPlot2dGpl);

	pMenuPlots->addSeparator();

	QAction *pSaveResults = new QAction("Save Results...", this);
	pSaveResults->setIcon(load_icon("res/icons/document-save-as.svg"));
	pMenuPlots->addAction(pSaveResults);

	// help menu
	QMenu *pMenuHelp = new QMenu("Help", this);

	QAction *pAbout = new QAction("About...", this);
	pAbout->setIcon(load_icon("res/icons/dialog-information.svg"));
	pMenuHelp->addAction(pAbout);


	m_pMenuBar->addMenu(pMenuFile);
	m_pMenuBar->addMenu(pMenuView);
	m_pMenuBar->addMenu(pMenuConvoActions);
	m_pMenuBar->addMenu(pMenuPlots);
	m_pMenuBar->addMenu(pMenuHelp);


	connect(pExit, &QAction::triggered, this, &ConvoDlg::accept);
	connect(pNew, &QAction::triggered, this, &ConvoDlg::New);
	connect(pLoad, &QAction::triggered, this, static_cast<void (ConvoDlg::*)()>(&ConvoDlg::Load));
	connect(pSave, &QAction::triggered, this, static_cast<void (ConvoDlg::*)()>(&ConvoDlg::Save));
	connect(pSaveAs, &QAction::triggered, this, &ConvoDlg::SaveAs);
	connect(pConvofit, &QAction::triggered, this, &ConvoDlg::SaveConvofit);
	connect(pActionLogY, &QAction::triggered, m_plotwrap.get(), &QwtPlotWrapper::ToggleLogY);
	connect(pActionStart, &QAction::triggered, this, &ConvoDlg::Start);
	connect(pActionStartFit, &QAction::triggered, this, &ConvoDlg::StartFit);
	connect(pActionDisp, &QAction::triggered, this, &ConvoDlg::StartDisp);
	connect(pExportPlot, &QAction::triggered, m_plotwrap.get(), &QwtPlotWrapper::SavePlot);
	connect(pExportPlot2d, &QAction::triggered, m_plotwrap2d.get(), &QwtPlotWrapper::SavePlot);
	connect(pExportPlotGpl, &QAction::triggered, m_plotwrap.get(), &QwtPlotWrapper::ExportGpl);
	connect(pExportPlot2dGpl, &QAction::triggered, m_plotwrap2d.get(), &QwtPlotWrapper::ExportGpl);
	connect(pSaveResults, &QAction::triggered, this, &ConvoDlg::SaveResult);
	connect(pAbout, &QAction::triggered, this, &ConvoDlg::ShowAboutDlg);

	this->layout()->setMenuBar(m_pMenuBar);
	// --------------------------------------------------------------------


	m_pSqwParamDlg = new SqwParamDlg(this, m_pSett);
	connect(this, &ConvoDlg::SqwLoaded, m_pSqwParamDlg, &SqwParamDlg::SqwLoaded);
	connect(m_pSqwParamDlg, &SqwParamDlg::SqwParamsChanged, this, &ConvoDlg::SqwParamsChanged);

	m_pFavDlg = new FavDlg(this, m_pSett);
	connect(m_pFavDlg, &FavDlg::ChangePos, this, &ConvoDlg::ChangePos);

	connect(btnBrowseCrys, &QToolButton::clicked, this, &ConvoDlg::browseCrysFiles);
	connect(btnBrowseRes, &QToolButton::clicked, this, &ConvoDlg::browseResoFiles);
	connect(btnBrowseSqw, &QPushButton::clicked, this, &ConvoDlg::browseSqwFiles);
	connect(btnBrowseScan, &QToolButton::clicked, this, &ConvoDlg::browseScanFiles);
	connect(btnBrowseAutosave, &QToolButton::clicked, this, &ConvoDlg::browseAutosaveFile);
	connect(btnFav, &QPushButton::clicked, this, &ConvoDlg::ShowFavourites);
	connect(btnSqwParams, &QToolButton::clicked, this, &ConvoDlg::showSqwParamDlg);

	connect(comboSqw, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), this, &ConvoDlg::SqwModelChanged);
	connect(editSqw, &QLineEdit::textChanged, this, &ConvoDlg::createSqwModel);
	connect(editScan, &QLineEdit::textChanged, this, &ConvoDlg::scanFileChanged);
	connect(editFilterCol, &QLineEdit::textChanged, this, [this]() -> void
	{ scanFileChanged(editScan->text()); });
	connect(editFilterVal, &QLineEdit::textChanged, this, [this]() -> void
	{ scanFileChanged(editScan->text()); });

	connect(editScale, &QLineEdit::textChanged, this, &ConvoDlg::scaleChanged);
	connect(editSlope, &QLineEdit::textChanged, this, &ConvoDlg::scaleChanged);
	connect(editOffs, &QLineEdit::textChanged, this, &ConvoDlg::scaleChanged);

	connect(btnStart, &QPushButton::clicked, this, &ConvoDlg::Start);
	connect(btnStartFit, &QPushButton::clicked, this, &ConvoDlg::StartFit);
	connect(btnStop, &QPushButton::clicked, this, &ConvoDlg::Stop);

	connect(checkScan, &QCheckBox::toggled, this, &ConvoDlg::scanCheckToggled);

	connect(pHK, &QAction::triggered, this, &ConvoDlg::ChangeHK);
	connect(pHL, &QAction::triggered, this, &ConvoDlg::ChangeHL);
	connect(pKL, &QAction::triggered, this, &ConvoDlg::ChangeKL);

	for(QDoubleSpinBox* pSpin : {spinStartH, spinStartK, spinStartL, spinStartE,
		spinStopH, spinStopK, spinStopL, spinStopE})
	{
		connect(pSpin, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			this, &ConvoDlg::UpdateCurFavPos);
	}

	LoadSettings();


#if BOOST_OS_MACOS
	// check if system python is available
	if(!tl::dir_exists("/Library/Frameworks/Python.framework")
		&& find_resource_dirs("Frameworks/Python.framework", false).size()==0)
	{
		QMessageBox::information(this, "Python Module",
			"The <i>Python</i> S(q,w) plugin module requires having the "
			"<a href=\"https://www.python.org/downloads/mac-osx/\">Python framework</a> installed.<br><br>"
			"Please also install the <i>numpy</i> and <i>scipy</i> packages using the following command:<br><br>"
			"<code>/Library/Frameworks/Python.framework/Versions/Current/bin/pip3 install numpy scipy</code>");
	}
#endif
}


ConvoDlg::~ConvoDlg()
{
	if(m_pth)
	{
		if(m_pth->joinable())
			m_pth->join();
		delete m_pth;
		m_pth = nullptr;
	}

	if(m_pSqwParamDlg) { delete m_pSqwParamDlg; m_pSqwParamDlg = nullptr; }
	if(m_pFavDlg) { delete m_pFavDlg; m_pFavDlg = nullptr; }
	if(m_pMenuBar) { delete m_pMenuBar; m_pMenuBar = nullptr; }

	if(m_pSqw) m_pSqw.reset();
	unload_sqw_plugins();
}


void ConvoDlg::SqwModelChanged(int)
{
	if(!m_bAllowSqwReinit) return;

	editSqw->clear();
	createSqwModel("");
}


void ConvoDlg::createSqwModel(const QString& qstrFile)
{
	if(!m_bAllowSqwReinit) return;

	if(m_pSqw)
	{
		m_pSqw.reset();
		emit SqwLoaded(std::vector<SqwBase::t_var>{}, nullptr);
	}

	std::string strSqwIdent = comboSqw->itemData(comboSqw->currentIndex()).toString().toStdString();
	std::string _strSqwFile = qstrFile.toStdString();
	tl::trim(_strSqwFile);
	std::string strSqwFile = find_file_in_global_paths(_strSqwFile);

	if(strSqwFile == "")
	{
		tl::log_warn("No S(q,w) config file given.");
		//return;
	}

	m_pSqw.reset();
	m_pSqw = construct_sqw(strSqwIdent, strSqwFile);
	if(!m_pSqw)
	{
		QMessageBox::critical(this, "Error", "Unknown S(q,w) model selected.");
		return;
	}

	if(m_pSqw && m_pSqw->IsOk())
	{
		emit SqwLoaded(m_pSqw->GetVars(), &m_pSqw->GetFitVars());
	}
	else
	{
		//QMessageBox::critical(this, "Error", "Could not create S(q,w).");
		tl::log_err("Could not create S(q,w).");
		return;
	}
}


void ConvoDlg::SqwParamsChanged(const std::vector<SqwBase::t_var>& vecVars,
	const std::vector<SqwBase::t_var_fit>* pvecVarsFit)
{
	if(!m_pSqw) return;

	m_pSqw->SetVars(vecVars);
	if(pvecVarsFit) m_pSqw->InitFitVars(*pvecVarsFit);

#ifndef NDEBUG
	// check: read parameters back in
	emit SqwLoaded(m_pSqw->GetVars(), &m_pSqw->GetFitVars());
#endif
}


/**
 * set a model parameter
 */
void ConvoDlg::SetSqwParam(const std::string& name, t_real_reso val)
{
	m_pSqw->SetVarIfAvail(name, tl::var_to_str(val));

	// read parameters back in to update paramters dialog
	emit SqwLoaded(m_pSqw->GetVars(), &m_pSqw->GetFitVars());
}


/**
 * set model parameters
 * [ ident, value, error ]
 */
void ConvoDlg::SetSqwParams(const std::vector<std::tuple<std::string, std::string, std::string>>& sqwparams)
{
	for(const auto& param : sqwparams)
	{
		m_pSqw->SetVarIfAvail(std::get<0>(param), std::get<1>(param));
		if(std::get<2>(param) != "")
			m_pSqw->SetErrIfAvail(std::get<0>(param), std::get<2>(param));
		//if(std::get<3>(param) != "")
		//	m_pSqw->SetRangeIfAvail(std::get<0>(param), std::get<3>(param));
	}

	/*auto vars = m_pSqw->GetVars();
	auto fitvars = m_pSqw->GetFitVars();
	for(const auto& var : vars)
		std::cout << std::get<0>(var) << " " << std::get<1>(var) << " " << std::get<2>(var) << std::endl;
	for(const auto& var : fitvars)
		std::cout << std::get<0>(var) << " " << std::get<1>(var) << " " << std::get<2>(var) << std::endl;*/

	// read parameters back in to update paramters dialog
	emit SqwLoaded(m_pSqw->GetVars(), &m_pSqw->GetFitVars());
}


/**
 * get the model parameters
 * [ ident, type, value, error, fit? ]
 */
ConvoDlg::t_sqwparams ConvoDlg::GetSqwParams(bool only_fitparams) const
{
	t_sqwparams params;

	std::vector<SqwBase::t_var> vars1 = m_pSqw->GetVars();
	std::vector<SqwBase::t_var_fit> vars2 = m_pSqw->GetFitVars();

	for(const SqwBase::t_var& var : vars1)
	{
		const std::string& strName = std::get<SQW_NAME>(var);
		std::string strErr, strRange;
		bool bFit = 0;

		// look for associated fit parameters: match with basic variable ident
		auto iterFit = std::find_if(vars2.begin(), vars2.end(),
			[&strName](const SqwBase::t_var_fit& varFit) -> bool
			{ return strName == std::get<0>(varFit); });

		if(iterFit != vars2.end())
		{
			strErr = std::get<1>(*iterFit);		// error
			strRange = std::get<3>(*iterFit);	// range
			bFit = std::get<2>(*iterFit);		// "is fit param" flag
		}

		if((only_fitparams && bFit) || !only_fitparams)
			params.emplace_back(std::make_tuple(strName, std::get<SQW_TYPE>(var), std::get<SQW_VAL>(var), strErr, bFit, strRange));
	}

	return params;
}

// -----------------------------------------------------------------------------


/**
 * clear plot curves
 */
void ConvoDlg::ClearPlot1D()
{
	static const std::vector<t_real> vecZero;
	if(!m_plotwrap) return;

	for(std::size_t iCurve=0; iCurve<CONVO_MAX_CURVES; ++iCurve)
		set_qwt_data<t_real_reso>()(*m_plotwrap, vecZero, vecZero, iCurve, false);
}


/**
 * start 1d or 2d convolutions
 */
void ConvoDlg::Start()
{
	if(check2dMap->isChecked())
		Start2D();
	else
		Start1D();
}


/**
 * stop running operations
 */
void ConvoDlg::Stop()
{
	m_atStop.store(true);
}


// -----------------------------------------------------------------------------


void ConvoDlg::ShowFavourites()
{
	focus_dlg(m_pFavDlg);
}

void ConvoDlg::UpdateCurFavPos()
{
	FavHklPos pos;
	pos.dhstart = spinStartH->value();
	pos.dkstart = spinStartK->value();
	pos.dlstart = spinStartL->value();
	pos.dEstart = spinStartE->value();
	pos.dhstop = spinStopH->value();
	pos.dkstop = spinStopK->value();
	pos.dlstop = spinStopL->value();
	pos.dEstop = spinStopE->value();

	m_pFavDlg->UpdateCurPos(pos);
}

void ConvoDlg::ChangePos(const struct FavHklPos& pos)
{
	spinStartH->setValue(pos.dhstart);
	spinStartK->setValue(pos.dkstart);
	spinStartL->setValue(pos.dlstart);
	spinStartE->setValue(pos.dEstart);
	spinStopH->setValue(pos.dhstop);
	spinStopK->setValue(pos.dkstop);
	spinStopL->setValue(pos.dlstop);
	spinStopE->setValue(pos.dEstop);
}


// -----------------------------------------------------------------------------


static void SwapSpin(QDoubleSpinBox* pS1, QDoubleSpinBox* pS2)
{
	double dVal = pS1->value();
	pS1->setValue(pS2->value());
	pS2->setValue(dVal);
}

void ConvoDlg::ChangeHK()
{
	SwapSpin(spinStartH, spinStartK);
	SwapSpin(spinStopH, spinStopK);
}
void ConvoDlg::ChangeHL()
{
	SwapSpin(spinStartH, spinStartL);
	SwapSpin(spinStopH, spinStopL);
}
void ConvoDlg::ChangeKL()
{
	SwapSpin(spinStartK, spinStartL);
	SwapSpin(spinStopK, spinStopL);
}


// -----------------------------------------------------------------------------

void ConvoDlg::scanCheckToggled(bool bChecked)
{
	if(bChecked)
		scanFileChanged(editScan->text());
}


void ConvoDlg::scanFileChanged(const QString& qstrFile)
{
	m_bUseScan = 0;
	if(!checkScan->isChecked())
		return;

	Filter filter;
	if(editFilterCol->text() != "")
	{
		filter.colEquals = std::make_pair(
			editFilterCol->text().toStdString(),
			editFilterVal->text().toStdString());
	}

	if(!load_scan_file(qstrFile.toStdString(), m_scan, checkFlip->isChecked(), filter))
	{
		tl::log_err("Cannot load scan(s).");
		return;
	}

	if(!m_scan.vecPoints.size())
	{
		tl::log_err("No points in scan(s).");
		return;
	}

	const ScanPoint& ptBegin = *m_scan.vecPoints.cbegin();
	const ScanPoint& ptEnd = *m_scan.vecPoints.crbegin();

	comboFixedK->setCurrentIndex(m_scan.bKiFixed ? 0 : 1);
	spinKfix->setValue(m_scan.dKFix);

	spinStartH->setValue(ptBegin.h);
	spinStartK->setValue(ptBegin.k);
	spinStartL->setValue(ptBegin.l);
	spinStartE->setValue(ptBegin.E / tl::get_one_meV<t_real>());

	spinStopH->setValue(ptEnd.h);
	spinStopK->setValue(ptEnd.k);
	spinStopL->setValue(ptEnd.l);
	spinStopE->setValue(ptEnd.E / tl::get_one_meV<t_real>());

	m_bUseScan = 1;
}


void ConvoDlg::scaleChanged()
{
	// get scan x axis
	bool bScanAxisFound = 0;
	int iScanAxisIdx = 0;
	std::string strScanVar = "";
	std::vector<std::vector<t_real>> vecAxes;
	std::tie(bScanAxisFound, iScanAxisIdx, strScanVar, vecAxes) = get_scan_axis<t_real>(
		true, comboAxis->currentIndex(), spinStepCnt->value(), EPS_RLU,
		spinStartH->value(), spinStopH->value(), spinStartK->value(), spinStopK->value(),
		spinStartL->value(), spinStopL->value(), spinStartE->value(), spinStopE->value());
	if(!bScanAxisFound)
	{
		tl::log_err("No scan variable found.");
		return;
	}


	t_real dScale = tl::str_to_var<t_real>(editScale->text().toStdString());
	t_real dSlope = tl::str_to_var<t_real>(editSlope->text().toStdString());
	t_real dOffs = tl::str_to_var<t_real>(editOffs->text().toStdString());

	m_vecScaledS.resize(m_vecS.size());
	for(std::size_t i=0; i<m_vecS.size(); ++i)
	{
		const t_real dXVal = vecAxes[iScanAxisIdx][i];
		m_vecScaledS[i] = dScale*(m_vecS[i] + dSlope*dXVal) + dOffs;
		if(m_vecScaledS[i] < 0.)
			m_vecScaledS[i] = 0.;
	}

	set_qwt_data<t_real_reso>()(*m_plotwrap, m_vecQ, m_vecScaledS, 0, false);
	set_qwt_data<t_real_reso>()(*m_plotwrap, m_vecQ, m_vecScaledS, 1, false);

	m_plotwrap->GetPlot()->replot();
}


// -----------------------------------------------------------------------------


void ConvoDlg::showSqwParamDlg()
{
	focus_dlg(m_pSqwParamDlg);
}


#include "libs/version.h"

void ConvoDlg::ShowAboutDlg()
{
	std::ostringstream ostrAbout;
	ostrAbout << "Takin/Convo version " << TAKIN_VER << ".\n";
	ostrAbout << "Written by Tobias Weber <tweber@ill.fr>,\n";
	ostrAbout << "2015 - 2021.\n";
	ostrAbout << "\n" << TAKIN_LICENSE("Takin/Convo");

	QMessageBox::about(this, "About Convo", ostrAbout.str().c_str());
}


#include "moc_ConvoDlg.cpp"
