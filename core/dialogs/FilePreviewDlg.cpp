/**
 * File Preview Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2015
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

#include "FilePreviewDlg.h"
#include <QGridLayout>
#include <QLabel>
#include <QPen>
#include "tlibs/file/loadinstr.h"
#include <iostream>
#include <memory>

using t_real = t_real_glob;


FilePreviewDlg::FilePreviewDlg(QWidget* pParent, const char* pcTitle, QSettings* pSett)
	: QFileDialog(pParent, pcTitle), m_pSettings(pSett)
{
	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);
	}

	m_pPlot.reset(new QwtPlot());
	m_plotwrap.reset(new QwtPlotWrapper(m_pPlot.get(), 2, true));

	QPen penCurve;
	penCurve.setColor(QColor(0,0,0x99));
	penCurve.setWidth(2);
	m_plotwrap->GetCurve(0)->setPen(penCurve);
	m_plotwrap->GetCurve(0)->setStyle(QwtPlotCurve::CurveStyle::Lines);
	m_plotwrap->GetCurve(0)->setTitle("Scan Curve");

	QPen penPoints;
	penPoints.setColor(QColor(0xff,0,0));
	penPoints.setWidth(4);
	m_plotwrap->GetCurve(1)->setPen(penPoints);
	m_plotwrap->GetCurve(1)->setStyle(QwtPlotCurve::CurveStyle::Dots);
	m_plotwrap->GetCurve(1)->setTitle("Scan Points");

	QSizePolicy spol(QSizePolicy::Preferred, QSizePolicy::Preferred);
	spol.setHorizontalStretch(1);
	spol.setVerticalStretch(1);
	m_pPlot->setSizePolicy(spol);

	// depends on qt/src/gui/dialogs/qfiledialog.ui using QGridLayout
	setOptions(QFileDialog::DontUseNativeDialog);
	QGridLayout *pLayout = (QGridLayout*)layout();
	if(pLayout)
	{
		pLayout->addWidget(m_pPlot.get(), pLayout->rowCount(), 0, 1, pLayout->columnCount());
		resize(size().width(), size().height()*1.25);
	}

	QObject::connect(this, &FilePreviewDlg::currentChanged,
		this, &FilePreviewDlg::FileSelected);
}

FilePreviewDlg::~FilePreviewDlg() {}

void FilePreviewDlg::ClearPlot()
{
	m_vecCts.clear();
	m_vecScn.clear();

	set_qwt_data<t_real>()(*m_plotwrap, m_vecScn, m_vecCts, 0, 0);
	set_qwt_data<t_real>()(*m_plotwrap, m_vecScn, m_vecCts, 1, 0);

	m_pPlot->setAxisTitle(QwtPlot::xBottom, "");
	m_pPlot->setAxisTitle(QwtPlot::yLeft, "");

	m_pPlot->replot();
}

void FilePreviewDlg::FileSelected(const QString& qstrFile)
{
	ClearPlot();

	tl::FileInstrBase<t_real> *pInstr = tl::FileInstrBase<t_real>::LoadInstr(qstrFile.toStdString().c_str());
	if(!pInstr) return;

	std::unique_ptr<tl::FileInstrBase<t_real>> _ptrInstr(pInstr);

	std::vector<std::string> vecScanVars = _ptrInstr->GetScannedVars();
	if(vecScanVars.size() == 0) return;
	m_vecCts = _ptrInstr->GetCol(_ptrInstr->GetCountVar());
	m_vecScn = _ptrInstr->GetCol(vecScanVars[0]);

	//std::copy(m_vecScn.begin(), m_vecScn.end(), std::ostream_iterator<t_real>(std::cout, " "));
	//std::cout << std::endl;

	set_qwt_data<t_real>()(*m_plotwrap, m_vecScn, m_vecCts, 0, 0);
	set_qwt_data<t_real>()(*m_plotwrap, m_vecScn, m_vecCts, 1, 0);

	m_pPlot->setAxisTitle(QwtPlot::xBottom, vecScanVars[0].c_str());
	m_pPlot->setAxisTitle(QwtPlot::yLeft, "Counts");

	m_pPlot->replot();
}


#include "moc_FilePreviewDlg.cpp"
