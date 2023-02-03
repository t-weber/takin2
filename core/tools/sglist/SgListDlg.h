/**
 * Space Group List Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date oct-2015
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

#ifndef __TAZ_SGLISTDLG_H__
#define __TAZ_SGLISTDLG_H__

#include <QDialog>
#include <QSettings>

#include "ui/ui_sglist.h"


class SgListDlg : public QDialog, Ui::SgListDlg
{ Q_OBJECT
protected:
	QSettings m_settings;

protected:
	void SetupSpacegroups();
	virtual void closeEvent(QCloseEvent* pEvt) override;

protected slots:
	void UpdateSG();
	void SGSelected(QListWidgetItem *pItem, QListWidgetItem *pItemPrev);
	void RecalcBragg();
	void SearchSG(const QString& qstr);

	void CalcTrafo();

public:
	SgListDlg(QWidget* pParent = nullptr);
	virtual ~SgListDlg();
};


#endif
