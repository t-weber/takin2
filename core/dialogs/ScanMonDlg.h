/**
 * Scan Monitor
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date Jul-2016
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

#ifndef __SCANMON_DLG_H__
#define __SCANMON_DLG_H__

#include <QDialog>
#include <QSettings>
#include "ui/ui_scanmon.h"

#include "libs/globals.h"
#include "libs/qt/qthelper.h"
#include "libs/qt/qwthelper.h"
#include "NetCacheDlg.h"


class ScanMonDlg : public QDialog, Ui::ScanMonDlg
{ Q_OBJECT
	protected:
		QSettings *m_pSettings = 0;

		t_real_glob m_dTimer = 0,
			m_dPreset = 0,
			m_dCounter = 0;

		std::vector<t_real_glob> m_vecX, m_vecY;
		std::unique_ptr<QwtPlotWrapper> m_plotwrap;

	protected:
		virtual void accept() override;

	public:
		ScanMonDlg(QWidget* pParent=0, QSettings *pSett=0);
		virtual ~ScanMonDlg() = default;

		void ClearPlot();
		void UpdatePlot(const std::string& strVals);

	public slots:
		void UpdateValue(const std::string& strKey, const CacheVal& val);
};

#endif
