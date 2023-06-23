/**
 * Real Space Parameters
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2014-2016
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __REAL_PARAMS_H__
#define __REAL_PARAMS_H__

#include <QDialog>
#include <QSettings>
#include <vector>

#include "ui/ui_real_params.h"
#include "libs/globals.h"
#include "libs/spacegroups/spacegroup.h"
#include "libs/spacegroups/latticehelper.h"
#include "tlibs/phys/lattice.h"
#include "AtomsDlg.h"


struct RealParams
{
	t_real_glob dMonoTT, dSampleTT, dAnaTT;
	t_real_glob dMonoT, dSampleT, dAnaT;
	t_real_glob dLenMonoSample, dLenSampleAna, dLenAnaDet;
};


class RealParamDlg : public QDialog, Ui::RealParamDlg
{ Q_OBJECT
	protected:
		QSettings *m_pSettings = 0;

		// lattice
		const tl::Lattice<t_real_glob>* m_pLatt = nullptr;

		// metric
		ublas::matrix<t_real_glob> m_matGCov, m_matGCont;

		bool m_bSamplePosSense = true;

	public:
		RealParamDlg(QWidget* pParent=0, QSettings* pSett=0);
		virtual ~RealParamDlg();

	public slots:
		void paramsChanged(const RealParams& parms);

		void CrystalChanged(const xtl::LatticeCommon<t_real_glob>&);
		void CalcVecs();
		void CalcCrystalRot();

		void SetSampleSense(bool bPos);

	protected:
		virtual void closeEvent(QCloseEvent *pEvt) override;
		virtual void showEvent(QShowEvent *pEvt) override;
		virtual void accept() override;
};

#endif
