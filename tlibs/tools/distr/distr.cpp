/**
 * random distributions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2016
 * @license GPLv2 or GPLv3
 */

#include "distr.h"
#include "string/string.h"
#include "math/math.h"
#include "log/log.h"

#include <sstream>
#include <qwt_picker_machine.h>

static const unsigned NUM_PTS = 512;


DistrDlg::DistrDlg(QWidget* pParent)
	: QDialog(pParent), m_settings("tobis_stuff", "distr_tool")
{
	setupUi(this);
	if(m_settings.contains("distr_dlg/geo"))
		restoreGeometry(m_settings.value("distr_dlg/geo").toByteArray());

	tl::DistrType tyDistrs[] = { tl::DistrType::NORMAL, tl::DistrType::LOGNORMAL,
		tl::DistrType::CAUCHY, tl::DistrType::CHI2,
		tl::DistrType::STUDENT, tl::DistrType::FISHER,
		tl::DistrType::EXP, tl::DistrType::BETA,
		tl::DistrType::GAMMA, tl::DistrType::LOGISTIC,
		tl::DistrType::BINOMIAL, tl::DistrType::HYPERGEOMETRIC, tl::DistrType::POISSON, };
	const char* pcDistrNames[] = { "Normal", "Log-Normal",
		"Cauchy", "Chi-Squared",
		"Student-t", "Fisher-f",
		"Exponential", "Beta",
		"Gamma", "Logistic",
		"Binomial", "Hypergeometric", "Poisson", };

	for(std::size_t iDistr=0; iDistr<sizeof(tyDistrs)/sizeof(tyDistrs[0]); ++iDistr)
	{
		QListWidgetItem *pItemDistr = new QListWidgetItem(pcDistrNames[iDistr]);
		pItemDistr->setData(Qt::UserRole, unsigned(tyDistrs[iDistr]));
		listDistr->addItem(pItemDistr);
	}

	QObject::connect(buttonBox, &QDialogButtonBox::accepted, this, &DistrDlg::accept);
	QObject::connect(listDistr, &QListWidget::currentItemChanged, this, &DistrDlg::DistrSelected);
	for(QLineEdit* pEdit : { editMinX, editMaxX, editParam1, editParam2, editParam3, editParam4 })
		QObject::connect(pEdit, &QLineEdit::textChanged, this, &DistrDlg::Calc);

	QPen penPDF, penCDF;
	penPDF.setColor(QColor(0,0,0xff));
	penCDF.setColor(QColor(0xff,0,0));
	penPDF.setWidth(2);
	penCDF.setWidth(2);

	m_pCurvePDF = new QwtPlotCurve("probability density");
	m_pCurveCDF = new QwtPlotCurve("cumulative distribution");
	m_pCurvePDF->setPen(penPDF);
	m_pCurveCDF->setPen(penCDF);
	m_pCurvePDF->setRenderHint(QwtPlotItem::RenderAntialiased, true);
	m_pCurveCDF->setRenderHint(QwtPlotItem::RenderAntialiased, true);
	m_pCurvePDF->attach(plot);
	m_pCurveCDF->attach(plot);

	m_pLegend = new QwtLegend();
	plot->insertLegend(m_pLegend, QwtPlot::BottomLegend);

	m_pPicker = new QwtPlotPicker(plot->xBottom, plot->yLeft, QwtPlotPicker::NoRubberBand,
		QwtPicker::AlwaysOff, plot->canvas());
	m_pPicker->setStateMachine(new QwtPickerTrackerMachine());
	QObject::connect(m_pPicker, &QwtPlotPicker::moved, this, &DistrDlg::PlotPicked);
	m_pPicker->setEnabled(1);

	listDistr->setCurrentItem(listDistr->item(0));
}

DistrDlg::~DistrDlg()
{}

void DistrDlg::DistrSelected(QListWidgetItem* pItem, QListWidgetItem*)
{
	if(!pItem) return;

	m_distrtype = tl::DistrType(pItem->data(Qt::UserRole).toUInt());
	Calc();
}

void DistrDlg::Calc()
{
	QLabel *pLabels[] = { labelParam1, labelParam2, labelParam3, labelParam4 };
	QLineEdit *pEdits[] = { editParam1, editParam2, editParam3, editParam4 };

	const t_real dMinX = tl::str_to_var<t_real>(editMinX->text().toStdString());
	const t_real dMaxX = tl::str_to_var<t_real>(editMaxX->text().toStdString());

	const t_real dParam1 = tl::str_to_var<t_real>(pEdits[0]->text().toStdString());
	const t_real dParam2 = tl::str_to_var<t_real>(pEdits[1]->text().toStdString());
	const t_real dParam3 = tl::str_to_var<t_real>(pEdits[2]->text().toStdString());
	const t_real dParam4 = tl::str_to_var<t_real>(pEdits[3]->text().toStdString());

	std::vector<const char*> vecParamNames;
	constexpr const char* strParam[] = { "Param 1", "Param 2", "Param 3", "Param 4" };

	try
	{
		if(m_distrtype == tl::DistrType::NORMAL)
		{
			using t_distr = tl::t_normal_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1, dParam2));
		}
		else if(m_distrtype == tl::DistrType::LOGNORMAL)
		{
			using t_distr = tl::t_lognormal_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1, dParam2));
		}
		else if(m_distrtype == tl::DistrType::CAUCHY)
		{
			using t_distr = tl::t_cauchy_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1, dParam2));
		}
		else if(m_distrtype == tl::DistrType::POISSON)
		{
			using t_distr = tl::t_poisson_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1));
		}
		else if(m_distrtype == tl::DistrType::BINOMIAL)
		{
			using t_distr = tl::t_binomial_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1, dParam2));
		}
		else if(m_distrtype == tl::DistrType::HYPERGEOMETRIC)
		{
			using t_distr = tl::t_hypergeo_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1, dParam2, dParam3));
		}
		else if(m_distrtype == tl::DistrType::CHI2)
		{
			using t_distr = tl::t_chi2_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1));
		}
		else if(m_distrtype == tl::DistrType::STUDENT)
		{
			using t_distr = tl::t_student_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1));
		}
		else if(m_distrtype == tl::DistrType::FISHER)
		{
			using t_distr = tl::t_fisher_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1, dParam2));
		}
		else if(m_distrtype == tl::DistrType::EXP)
		{
			using t_distr = tl::t_exp_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1));
		}
		else if(m_distrtype == tl::DistrType::BETA)
		{
			using t_distr = tl::t_beta_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1, dParam2));
		}
		else if(m_distrtype == tl::DistrType::GAMMA)
		{
			using t_distr = tl::t_gamma_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1, dParam2));
		}
		else if(m_distrtype == tl::DistrType::LOGISTIC)
		{
			using t_distr = tl::t_logistic_dist<t_real>;
			vecParamNames = t_distr::GetParamNamesVec();
			m_pDistr.reset(new t_distr(dParam1, dParam2));
		}
		else
		{
			return;
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}

	try
	{
		std::size_t iNumParams = vecParamNames.size();
		std::size_t iEdit = 0;
		for(iEdit=0; iEdit<iNumParams; ++iEdit)
		{
			pLabels[iEdit]->setEnabled(1);
			pEdits[iEdit]->setEnabled(1);
			pLabels[iEdit]->setText(vecParamNames[iEdit]);
		}
		for(; iEdit<sizeof(pEdits)/sizeof(pEdits[0]); ++iEdit)
		{
			pLabels[iEdit]->setEnabled(0);
			pEdits[iEdit]->setEnabled(0);
			pLabels[iEdit]->setText(strParam[iEdit]);
		}

		m_vecX.clear();
		m_vecPDF.clear();
		m_vecCDF.clear();

		if(!!m_pDistr)
		{
			m_vecX = tl::linspace<t_real>(dMinX, dMaxX, NUM_PTS);
			m_vecPDF.resize(m_vecX.size());
			m_vecCDF.resize(m_vecX.size());
			for(std::size_t iX=0; iX<m_vecX.size(); ++iX)
			{
				m_vecPDF[iX] = m_pDistr->pdf(m_vecX[iX]);
				m_vecCDF[iX] = m_pDistr->cdf(m_vecX[iX]);
			}
		}

		m_pCurvePDF->setRawSamples(m_vecX.data(), m_vecPDF.data(), m_vecX.size());
		m_pCurveCDF->setRawSamples(m_vecX.data(), m_vecCDF.data(), m_vecX.size());
		plot->replot();
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}

void DistrDlg::PlotPicked(const QPointF& pt)
{
	t_real dX = t_real(pt.x());

	std::ostringstream ostr;
	ostr << "x = " << dX;

	try
	{
		if(!!m_pDistr)
		{
			t_real dPDF = m_pDistr->pdf(dX);
			t_real dCDF = m_pDistr->cdf(dX);

			ostr << ", pdf(x) = " << dPDF;
			ostr << ", cdf(x) = " << dCDF;
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}

	labelStatus->setText(ostr.str().c_str());
}

void DistrDlg::accept()
{
	m_settings.setValue("distr_dlg/geo", saveGeometry());
	QDialog::accept();
}
