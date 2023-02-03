/**
 * Spurion Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 26-may-2014
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

#ifndef __SPURION_DLG_H__
#define __SPURION_DLG_H__

#include <QDialog>
#include <QSettings>
#include <vector>
#include <memory>
#include "ui/ui_spurions.h"
#include "RecipParamDlg.h"
#include "libs/qt/qthelper.h"
#include "libs/qt/qwthelper.h"
#include "libs/globals.h"


class SpurionDlg : public QDialog, Ui::SpurionDlg
{ Q_OBJECT
	protected:
		QSettings *m_pSettings = 0;
		t_real_glob m_dEi=0., m_dEf=0.;

		std::vector<t_real_glob> m_vecQ, m_vecE;
		std::unique_ptr<QwtPlotWrapper> m_plotwrap;

	public:
		SpurionDlg(QWidget* pParent=0, QSettings *pSett=0);
		virtual ~SpurionDlg();

	protected slots:
		void ChangedKiKfMode();
		void Calc();

		void CalcInel();
		void CalcBragg();

		void cursorMoved(const QPointF& pt);
		void SaveTable();

	public slots:
		void paramsChanged(const RecipParams& parms);

	protected:
		virtual void showEvent(QShowEvent *pEvt) override;
		virtual void accept() override;
};


#endif
