/**
 * 3D Ellipsoid Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date may-2013, 29-apr-2014
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

#ifndef __RESO_ELLI_DLG_3D__
#define __RESO_ELLI_DLG_3D__

#include <QDialog>
#include <QSettings>
#include <QComboBox>

#include <vector>

#include "EllipseDlg.h"
#include "libs/plotgl.h"
#include "tlibs/math/linalg.h"
#include "tools/res/ellipse.h"
#include "tools/res/defs.h"


class EllipseDlg3D : public QDialog
{Q_OBJECT
	protected:
		std::vector<PlotGl*> m_pPlots;
		std::vector<Ellipsoid3d<t_real_reso>> m_elliProj;
		std::vector<Ellipsoid3d<t_real_reso>> m_elliSlice;

		QComboBox *m_pComboCoord = nullptr;
		QSettings *m_pSettings = nullptr;

		ublas::matrix<t_real_reso> m_reso, m_resoHKL, m_resoOrient;
		ublas::vector<t_real_reso> m_reso_v = ublas::zero_vector<t_real_reso>(4),
			m_reso_vHKL = ublas::zero_vector<t_real_reso>(4),
			m_reso_vOrient = ublas::zero_vector<t_real_reso>(4);
		t_real_reso m_reso_s = 0;
		ublas::vector<t_real_reso> m_Q_avg, m_Q_avgHKL, m_Q_avgOrient;
		ResoAlgo m_algo = ResoAlgo::UNKNOWN;

	protected:
		ublas::vector<t_real_reso>
		ProjRotatedVec(const ublas::matrix<t_real_reso>& rot,
			const ublas::vector<t_real_reso>& vec);

	public:
		EllipseDlg3D(QWidget* pParent, QSettings *pSett=0);
		virtual ~EllipseDlg3D();

	protected:
		virtual void hideEvent(QHideEvent*) override;
		virtual void showEvent(QShowEvent*) override;
		virtual void closeEvent(QCloseEvent*) override;
		virtual void keyPressEvent(QKeyEvent*) override;
		virtual void accept() override;

	public slots:
		void SetParams(const EllipseDlgParams& params);
		void Calc();
};

#endif
