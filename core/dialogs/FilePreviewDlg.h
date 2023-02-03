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

#ifndef __FILE_PREV_DLG__
#define __FILE_PREV_DLG__

#include <vector>
#include <memory>

#include <QFileDialog>
#include <QString>
#include <QSettings>

#include <qwt_plot.h>
#include <qwt_plot_curve.h>

#include "libs/qt/qthelper.h"
#include "libs/qt/qwthelper.h"
#include "libs/globals.h"


class FilePreviewDlg : public QFileDialog
{ Q_OBJECT
	protected:
		QSettings *m_pSettings = nullptr;

		std::unique_ptr<QwtPlot> m_pPlot;
		std::unique_ptr<QwtPlotWrapper> m_plotwrap;
		std::vector<t_real_glob> m_vecScn, m_vecCts;

	protected:
		void ClearPlot();

	public:
		FilePreviewDlg(QWidget* pParent, const char* pcTitle, QSettings* pSett=0);
		virtual ~FilePreviewDlg();

	protected slots:
		void FileSelected(const QString& strFile);
};

#endif
