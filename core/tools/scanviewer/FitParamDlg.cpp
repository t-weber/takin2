/**
 * Fit Parameters Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2016
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

#include "FitParamDlg.h"

FitParamDlg::FitParamDlg(QWidget* pParent, QSettings *pSett)
	: QDialog(pParent), m_pSettings(pSett)
{
	this->setupUi(this);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);

		if(m_pSettings->contains("fit_params/geo"))
			restoreGeometry(m_pSettings->value("fit_params/geo").toByteArray());
	}
}

FitParamDlg::~FitParamDlg()
{}


void FitParamDlg::SetBold(QLabel* pLab, bool bBold)
{
	QFont fontLabel = pLab->font();
	fontLabel.setBold(bBold);
	pLab->setFont(fontLabel);
}

void FitParamDlg::UnsetAllBold()
{
	for(QLabel* pLab : {labelAmp, labelSig, labelHWHM, labelX0, labelOffs, labelSlope, labelFreq, labelPhase})
		SetBold(pLab, 0);
}


void FitParamDlg::accept()
{
	if(m_pSettings)
		m_pSettings->setValue("fit_params/geo", saveGeometry());

	QDialog::accept();
}


#include "moc_FitParamDlg.cpp"
