/**
 * random distributions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2016
 * @license GPLv2 or GPLv3
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
