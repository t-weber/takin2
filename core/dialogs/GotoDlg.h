/**
 * Goto Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 15-oct-2014
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

#ifndef __GOTO_DLG_H__
#define __GOTO_DLG_H__

#include <QDialog>
#include <QSettings>
#include "ui/ui_goto.h"

#include "tlibs/phys/lattice.h"
#include "tlibs/math/linalg.h"
#include "tlibs/file/prop.h"
#include "libs/globals.h"
#include "libs/globals_qt.h"
#include "tools/taz/tasoptions.h"
#include "RecipParamDlg.h"


class GotoDlg : public QDialog, Ui::GotoDlg
{ Q_OBJECT
	private:
		bool m_bMonoAnaOk = 0;
		bool m_bSampleOk = 0;

	public:
		GotoDlg(QWidget* pParent=0, QSettings* pSett=0);
		virtual ~GotoDlg();

	protected:
		QSettings *m_pSettings = 0;

		ublas::vector<t_real_glob> m_vec1, m_vec2;
		tl::Lattice<t_real_glob> m_lattice;
		RecipParams m_paramsRecip;
		bool m_bHasParamsRecip = 0;

		t_real_glob m_dMono2Theta;
		t_real_glob m_dSample2Theta;
		t_real_glob m_dSampleTheta;
		t_real_glob m_dAna2Theta;

		t_real_glob m_dAna = 3.355;
		t_real_glob m_dMono = 3.355;

		bool m_bSenseM=0, m_bSenseS=1, m_bSenseA=0;

	public:
		void ClearList();

	protected:
		bool GotoPos(QListWidgetItem* pItem, bool bApply);
		bool ApplyCurPos();

	protected slots:
		void EditedKiKf();
		void EditedE();
		void EditedAngles();
		void GetCurPos();

		// list
		void AddPosToList();
		void RemPosFromList();
		void LoadList();
		void SaveList();
		void ListItemSelected();
		void ListItemDoubleClicked(QListWidgetItem*);

	public slots:
		void CalcMonoAna();
		void CalcSample();

		void AddPosToList(t_real_glob dh, t_real_glob dk, t_real_glob dl, t_real_glob dki, t_real_glob dkf);

	public:
		void SetLattice(const tl::Lattice<t_real_glob>& lattice) { m_lattice = lattice; }
		void SetScatteringPlane(const ublas::vector<t_real_glob>& vec1, const ublas::vector<t_real_glob>& vec2)
		{ m_vec1 = vec1; m_vec2 = vec2; }

		void SetD(t_real_glob dMono, t_real_glob dAna) { m_dMono = dMono; m_dAna = dAna; }
		void SetSenses(bool bM, bool bS, bool bA)
		{ m_bSenseM=bM; m_bSenseS=bS; m_bSenseA=bA; }

		void SetMonoSense(bool bSense) { m_bSenseM = bSense; }
		void SetSampleSense(bool bSense) { m_bSenseS = bSense; }
		void SetAnaSense(bool bSense) { m_bSenseA = bSense; }

		bool GotoPos(unsigned int iItem);

	protected slots:
		void ButtonBoxClicked(QAbstractButton* pBtn);
	public slots:
		void RecipParamsChanged(const RecipParams&);

	signals:
		void vars_changed(const CrystalOptions& crys, const TriangleOptions& triag);

	protected:
		virtual void showEvent(QShowEvent *pEvt) override;

	public:
		void Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot);
		void Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot);
};

#endif
