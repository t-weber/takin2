/**
 * Neutron Properties Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jul-2013, 28-may-2014
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

#ifndef __NEUTRON_DLG_H__
#define __NEUTRON_DLG_H__

#include <QDialog>
#include <QSettings>
#include "ui/ui_neutrons.h"

#include "libs/globals.h"
#include "RecipParamDlg.h"


class NeutronDlg : public QDialog, Ui::NeutronDlg
{ Q_OBJECT
	protected:
		QSettings *m_pSettings = 0;

		t_real_glob m_dExtKi = 0.;
		t_real_glob m_dExtKf = 0.;

	protected:
		void setupConstants();

	public:
		NeutronDlg(QWidget* pParent=0, QSettings* pSett=0);
		virtual ~NeutronDlg();

	protected slots:
		void CalcNeutronLam();
		void CalcNeutronk();
		void CalcNeutronv();
		void CalcNeutronE();
		void CalcNeutronOm();
		void CalcNeutronF();
		void CalcNeutronT();
		void CalcNeutronTau();

		void RecipThetaEdited();
		void RecipTwoThetaEdited();
		void RecipKEdited();
		void RecipLamEdited();
		void RecipGEdited();
		void RecipDEdited();

		void CalcBraggRecip();

		void EnableRecipEdits();

		void Eval(const QString&);

		void SetExtKi();
		void SetExtKf();

	public slots:
		void paramsChanged(const RecipParams& parms);

	protected:
		virtual void showEvent(QShowEvent *pEvt) override;
		virtual void accept() override;

		static void SetEditTT(QLineEdit *pEditT, QLineEdit *pEditTT);
		static void SetEditT(QLineEdit *pEditT, QLineEdit *pEditTT);

		static void SetEditK(QLineEdit *pEditLam, QLineEdit *pEditK);
		static void SetEditLam(QLineEdit *pEditLam, QLineEdit *pEditK);

		static void SetEditG(QLineEdit *pEditG, QLineEdit *pEditD);
		static void SetEditD(QLineEdit *pEditG, QLineEdit *pEditD);
};


#endif
