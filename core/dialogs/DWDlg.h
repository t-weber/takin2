/**
 * Scattering factors dialog (e.g. Debye-Waller factor)
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013, jan-2015
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

#ifndef __DWDLG_H__
#define __DWDLG_H__

#include <QDialog>
#include <QSettings>
#include "ui/ui_dw.h"

#include <vector>
#include <memory>
#include "libs/qt/qthelper.h"
#include "libs/qt/qwthelper.h"
#include "libs/globals.h"
#include <qwt_legend.h>


class DWDlg : public QDialog, Ui::DWDlg
{ Q_OBJECT
protected:
	QSettings *m_pSettings = nullptr;

	// dw stuff
	std::vector<t_real_glob> m_vecQ, m_vecDeb;
	std::unique_ptr<QwtPlotWrapper> m_plotwrapDW;

	// ana stuff
	std::vector<t_real_glob> m_veckf, m_vecInt;
	std::unique_ptr<QwtPlotWrapper> m_plotwrapAna;

	// bose stuff
	std::vector<t_real_glob> m_vecBoseE, m_vecBoseIntPos, m_vecBoseIntNeg;
	std::unique_ptr<QwtPlotWrapper> m_plotwrapBose;
	QwtLegend *m_pLegendBose = nullptr;

	// lorentz stuff
	std::vector<t_real_glob> m_vecLor2th, m_vecLor;
	std::unique_ptr<QwtPlotWrapper> m_plotwrapLor;

protected:
	virtual void showEvent(QShowEvent *pEvt) override;
	virtual void accept() override;

protected slots:
	void cursorMoved(const QPointF& pt);
	void CalcDW();
	void CalcAna();
	void CalcBose();
	void CalcLorentz();

public:
	DWDlg(QWidget* pParent = nullptr, QSettings *pSettings = nullptr);
	virtual ~DWDlg();
};

#endif
