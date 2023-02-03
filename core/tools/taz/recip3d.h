/**
 * Scattering Triangle Tool
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date mar-2014
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

#ifndef __TAZ_RECIP_3D__
#define __TAZ_RECIP_3D__

#include <QDialog>
#include <QStatusBar>
#include <boost/optional.hpp>

#include "tlibs/math/linalg.h"
#include "tlibs/math/geo.h"
#include "tlibs/phys/lattice.h"

#include "libs/plotgl.h"
#include "libs/spacegroups/spacegroup.h"
#include "libs/spacegroups/latticehelper.h"
#include "libs/globals.h"

#include "dialogs/RecipParamDlg.h"


class Recip3DDlg : public QDialog
{Q_OBJECT
protected:
	QSettings* m_pSettings = nullptr;
	QStatusBar *m_pStatus = nullptr;
	PlotGl* m_pPlot = nullptr;
	t_real_glob m_dMaxPeaks = 4.;
	t_real_glob m_dPlaneDistTolerance = 0.01;

	boost::optional<std::size_t> m_iQIdx;
	ublas::vector<t_real_glob> m_vecQ, m_vecQ_rlu;

public:
	Recip3DDlg(QWidget* pParent, QSettings* = 0);
	virtual ~Recip3DDlg();

	void CalcPeaks(const xtl::LatticeCommon<t_real_glob>& recipcommon);

	void SetPlaneDistTolerance(t_real_glob dTol) { m_dPlaneDistTolerance = dTol; }
	void SetMaxPeaks(t_real_glob dMax) { m_dMaxPeaks = dMax; }

protected:
	virtual void hideEvent(QHideEvent*) override;
	virtual void showEvent(QShowEvent*) override;
	virtual void closeEvent(QCloseEvent*) override;
	virtual void keyPressEvent(QKeyEvent*) override;

public slots:
	void RecipParamsChanged(const RecipParams&);
};


#endif
