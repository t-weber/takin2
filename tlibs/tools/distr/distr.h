/**
 * random distributions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2016
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __DISTR_H__
#define __DISTR_H__

#include <QMainWindow>
#include <QSettings>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_legend.h>
#include <qwt_plot_picker.h>

#include <vector>

#include "math/distr.h"
#include "ui_distr.h"

typedef double t_real;

class DistrDlg : public QDialog, Ui::DistrDlg
{
	protected:
		QSettings m_settings;

		tl::DistrType m_distrtype = tl::DistrType::NONE;
		std::unique_ptr<tl::DistrBase<t_real>> m_pDistr;

		std::vector<t_real> m_vecX, m_vecPDF, m_vecCDF;
		QwtPlotCurve *m_pCurvePDF=nullptr, *m_pCurveCDF=nullptr;
		QwtLegend *m_pLegend=nullptr;
		QwtPlotPicker* m_pPicker=nullptr;

	public:
		DistrDlg(QWidget *pParent=nullptr);
		virtual ~DistrDlg();

	protected:
		void DistrSelected(QListWidgetItem*, QListWidgetItem*);
		void Calc();
		void PlotPicked(const QPointF& pt);

		virtual void accept() override;
};

#endif
