/**
 * Atom Positions Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
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

#ifndef __TAKIN_ATOMS_DLG_H__
#define __TAKIN_ATOMS_DLG_H__

#include <QDialog>
#include <QSettings>
#include "tlibs/helper/boost_hacks.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <string>

#include "libs/globals.h"
#include "libs/globals_qt.h"
#include "libs/spacegroups/latticehelper.h"

#include "ui/ui_atoms.h"
namespace ublas = boost::numeric::ublas;



class AtomsDlg : public QDialog, Ui::AtomsDlg
{ Q_OBJECT
protected:
	QSettings *m_pSettings = nullptr;
	bool m_bEnableSpin = 0;
	std::string m_strErr;

protected:
	virtual void closeEvent(QCloseEvent*) override;
	void SendApplyAtoms();
	void CheckAtoms();

protected slots:
	void ButtonBoxClicked(QAbstractButton* pBtn);
	void RemoveAtom();
	void AddAtom();
	void AtomCellChanged(int iRow, int iCol);

public:
	AtomsDlg(QWidget* pParent = nullptr, QSettings *pSettings = nullptr,
		bool bEnableSpin=0);
	virtual ~AtomsDlg();

	void SetAtoms(const std::vector<xtl::AtomPos<t_real_glob>>& vecAtoms);
	std::string GetErrorString() const { return m_strErr; }
	bool ShowPossibleErrorDlg();

signals:
	void ApplyAtoms(const std::vector<xtl::AtomPos<t_real_glob>>& vecAtoms);
};


#endif
