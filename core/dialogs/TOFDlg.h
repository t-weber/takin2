/**
 * Component calculations (formerly only TOF-specific)
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2017
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

#ifndef __TOFCALC_DLG_H__
#define __TOFCALC_DLG_H__

#include <QDialog>
#include <QSettings>
#include "ui/ui_tof.h"

#include "libs/globals.h"

#include <vector>
#include <string>


class TOFDlg : public QDialog, Ui::TofCalcDlg
{ Q_OBJECT
	protected:
		QSettings *m_pSettings = nullptr;

		std::vector<QRadioButton*> m_vecRadioBoxes;
		std::vector<std::string> m_vecRadioNames;

		std::vector<QCheckBox*> m_vecCheckBoxes;
		std::vector<std::string> m_vecCheckNames;

		std::vector<QLineEdit*> m_vecEditBoxes;
		std::vector<std::string> m_vecEditNames;

	public:
		TOFDlg(QWidget* pParent=0, QSettings* pSett=0);
		virtual ~TOFDlg() = default;

	protected:
		virtual void showEvent(QShowEvent *pEvt) override;
		virtual void accept() override;
		void ReadLastConfig();
		void WriteLastConfig();

	protected slots:
		void CalcChopper();
		void EnableChopperEdits();

		void CalcDiv();
		void EnableDivEdits();

		void CalcSel();
		void EnableSelEdits();

		void CalcFoc();
};


#endif
