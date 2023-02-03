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

#include "SgListDlg.h"
#include <sstream>
#include "tlibs/string/string.h"
#include "tlibs/math/linalg.h"
#include "libs/spacegroups/spacegroup.h"
#include "libs/globals.h"

using t_real = t_real_glob;
namespace ublas = tl::ublas;


SgListDlg::SgListDlg(QWidget *pParent)
	: QDialog(pParent, Qt::WindowTitleHint|Qt::WindowCloseButtonHint|Qt::WindowMinMaxButtonsHint),
		m_settings("takin", "sglist")
{
	setupUi(this);
	QFont font;
	if(m_settings.contains("main/font_gen") && font.fromString(m_settings.value("main/font_gen", "").toString()))
		setFont(font);


	SetupSpacegroups();

	QObject::connect(listSGs, &QListWidget::currentItemChanged, this, &SgListDlg::SGSelected);
	QObject::connect(checkMatrices, &QCheckBox::toggled, this, &SgListDlg::UpdateSG);

	for(QSpinBox* pSpin : {spinH, spinK, spinL})
		QObject::connect(pSpin, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
			this, &SgListDlg::RecalcBragg);

	for(QDoubleSpinBox* pSpin : {spinX, spinY, spinZ, spinW})
		QObject::connect(pSpin, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			this, &SgListDlg::CalcTrafo);

	QObject::connect(editFilter, &QLineEdit::textEdited, this, &SgListDlg::SearchSG);


	if(m_settings.contains("sglist/geo"))
		restoreGeometry(m_settings.value("sglist/geo").toByteArray());
}

SgListDlg::~SgListDlg()
{}


void SgListDlg::closeEvent(QCloseEvent* pEvt)
{
	m_settings.setValue("sglist/geo", saveGeometry());
}


static QListWidgetItem* create_header_item(const char *pcTitle)
{
	QListWidgetItem *pHeaderItem = new QListWidgetItem(pcTitle);
	pHeaderItem->setTextAlignment(Qt::AlignHCenter);

	QFont fontHeader = pHeaderItem->font();
	fontHeader.setBold(1);
	pHeaderItem->setFont(fontHeader);

	QBrush brushHeader = pHeaderItem->foreground();
	brushHeader.setColor(QColor(0xff, 0xff, 0xff));
	pHeaderItem->setForeground(brushHeader);
	pHeaderItem->setBackground(QColor(0x65, 0x65, 0x65));

	pHeaderItem->setData(Qt::UserRole, 1000);

	return pHeaderItem;
}

void SgListDlg::SetupSpacegroups()
{
	listSGs->clear();

	std::shared_ptr<const xtl::SpaceGroups<t_real>> sgs = xtl::SpaceGroups<t_real>::GetInstance();
	const xtl::SpaceGroups<t_real>::t_vecSpaceGroups* pvecSG = sgs->get_space_groups_vec();

	// prevents double insertion of headers if two space groups have the same number
	bool bAlreadySeen[7] = { 0, 0, 0, 0, 0, 0, 0 };
	const char** pcHeader = xtl::get_crystal_system_names(1);
	const unsigned int *piStartNr = xtl::get_crystal_system_start_indices();
	const QColor itemColsBackground[] = {QColor(0xff, 0xff, 0xff), QColor(0xee, 0xee, 0xee)};
	const QColor itemColsForeground[] = {QColor(0x00, 0x00, 0x00), QColor(0x00, 0x00, 0x00)};

	for(unsigned int iSG=0; iSG<pvecSG->size(); ++iSG)
	{
		const xtl::SpaceGroup<t_real>* psg = pvecSG->at(iSG);
		unsigned int iSgNr = psg->GetNr();

		// crystal system headers
		for(unsigned int iCrystSys=0; iCrystSys<7; ++iCrystSys)
		{
			if(iSgNr==piStartNr[iCrystSys] && !bAlreadySeen[iCrystSys])
			{
				listSGs->addItem(create_header_item(pcHeader[iCrystSys]));
				bAlreadySeen[iCrystSys] = 1;
				break;
			}
		}

		std::ostringstream ostrSg;
		ostrSg << "No. " << iSgNr << ": " << psg->GetName();

		QListWidgetItem* pItem = new QListWidgetItem(ostrSg.str().c_str());
		pItem->setData(Qt::UserRole, iSG);
		pItem->setBackground(itemColsBackground[iSgNr % (sizeof(itemColsBackground)/sizeof(itemColsBackground[0]))]);
		pItem->setForeground(itemColsForeground[iSgNr % (sizeof(itemColsForeground)/sizeof(itemColsForeground[0]))]);
		listSGs->addItem(pItem);
	}
}

void SgListDlg::UpdateSG()
{
	SGSelected(listSGs->currentItem(), nullptr);
}

void SgListDlg::SGSelected(QListWidgetItem *pItem, QListWidgetItem*)
{
	listSymOps->clear();
	listTrafo->clear();
	for(QLineEdit *pEdit : {editHM, /*editHall,*/ editLaue, editNr/*, editAxisSym*/})
		pEdit->setText("");
	if(!pItem) return;


	std::shared_ptr<const xtl::SpaceGroups<t_real>> sgs = xtl::SpaceGroups<t_real>::GetInstance();
	const xtl::SpaceGroups<t_real>::t_vecSpaceGroups* pvecSG = sgs->get_space_groups_vec();

	// header selected?
	unsigned int iSG = pItem->data(Qt::UserRole).toUInt();
	if(iSG >= pvecSG->size())
		return;

	const xtl::SpaceGroup<t_real>* psg = pvecSG->at(iSG);
	unsigned int iSgNr = psg->GetNr();

	const std::string& strHM = psg->GetName();
	const std::string& strPointGroup = psg->GetPointGroup();
	const std::string& strLaue = psg->GetLaueGroup();
	const std::string& strCrysSys = psg->GetCrystalSystemName();

	editNr->setText(tl::var_to_str(iSgNr).c_str());
	editHM->setText(strHM.c_str());
	//editHall->setText(psg.symbol_hall().c_str());

	std::string strPtGr = "PG: " + strPointGroup;
	if(strLaue != "")
		strPtGr += ", LG: " + strLaue;
	strPtGr += " (" + strCrysSys + ")";
	editLaue->setText(strPtGr.c_str());

	bool bShowMatrices = checkMatrices->isChecked();

	// all trafos
	const std::vector<xtl::SpaceGroup<t_real>::t_mat>& vecTrafos = psg->GetTrafos();
	{
		std::ostringstream ostr;
		ostr << "All Symmetry Operations (" << vecTrafos.size() << ")";
		listSymOps->addItem(create_header_item(ostr.str().c_str()));

		for(unsigned int iSymOp=0; iSymOp<vecTrafos.size(); ++iSymOp)
		{
			std::string strDesc;
			if(bShowMatrices)
				strDesc = xtl::print_matrix(vecTrafos[iSymOp]);
			else
				strDesc = xtl::get_trafo_desc(vecTrafos[iSymOp]);

	                QListWidgetItem* pItem = new QListWidgetItem(strDesc.c_str());
        	        pItem->setData(Qt::UserRole, iSymOp);
			listSymOps->addItem(pItem);
		}
	}


	// primitive trafos
	const std::vector<unsigned int>& vecPrim = psg->GetPrimTrafos();

	if(vecPrim.size())
	{
		std::ostringstream ostr;
		ostr << "Primitive Symmetry Operations (" << (vecPrim.size()) << ")";
		listSymOps->addItem(create_header_item(ostr.str().c_str()));
		for(unsigned int iSymOp=0; iSymOp<vecPrim.size(); ++iSymOp)
		{
			if(bShowMatrices)
				listSymOps->addItem(xtl::print_matrix(vecTrafos[vecPrim[iSymOp]]).c_str());
			else
				listSymOps->addItem(xtl::get_trafo_desc(vecTrafos[vecPrim[iSymOp]]).c_str());
		}
	}


	// inverting trafos
	const std::vector<unsigned int>& vecInv = psg->GetInvTrafos();

	if(vecInv.size())
	{
		std::ostringstream ostr;
		ostr << "Inverting Symmetry Operations (" << (vecInv.size()) << ")";
		listSymOps->addItem(create_header_item(ostr.str().c_str()));
		for(unsigned int iSymOp=0; iSymOp<vecInv.size(); ++iSymOp)
		{
			if(bShowMatrices)
				listSymOps->addItem(xtl::print_matrix(vecTrafos[vecInv[iSymOp]]).c_str());
			else
				listSymOps->addItem(xtl::get_trafo_desc(vecTrafos[vecInv[iSymOp]]).c_str());
		}
	}

	// centering trafos
	const std::vector<unsigned int>& vecCenter = psg->GetCenterTrafos();

	if(vecCenter.size())
	{
		std::ostringstream ostr;
		ostr << "Centering Symmetry Operations (" << (vecCenter.size()) << ")";
		listSymOps->addItem(create_header_item(ostr.str().c_str()));
		for(unsigned int iSymOp=0; iSymOp<vecCenter.size(); ++iSymOp)
		{
			if(bShowMatrices)
				listSymOps->addItem(xtl::print_matrix(vecTrafos[vecCenter[iSymOp]]).c_str());
			else
				listSymOps->addItem(xtl::get_trafo_desc(vecTrafos[vecCenter[iSymOp]]).c_str());
		}
	}

	const std::vector<unsigned int>& vecTrans = psg->GetTransTrafos();

	if(vecTrans.size())
	{
		std::ostringstream ostr;
		ostr << "Translating Symmetry Operations (" << (vecTrans.size()) << ")";
		listSymOps->addItem(create_header_item(ostr.str().c_str()));
		for(unsigned int iSymOp=0; iSymOp<vecTrans.size(); ++iSymOp)
		{
			if(bShowMatrices)
				listSymOps->addItem(xtl::print_matrix(vecTrafos[vecTrans[iSymOp]]).c_str());
			else
				listSymOps->addItem(xtl::get_trafo_desc(vecTrafos[vecTrans[iSymOp]]).c_str());
		}
	}

/*
	// screw axes (rotation around and translation along axis)
	// and glide planes (reflection at and translation parallel to plane)
	const std::vector<unsigned int>& vecScrews = psg->GetScrewsNGlides();

	if(vecScrews.size())
	{
		std::ostringstream ostr;
		ostr << "Screw/Glide Symmetry Operations (" << (vecScrews.size()) << ")";
		listSymOps->addItem(create_header_item(ostr.str().c_str()));
		for(unsigned int iSymOp=0; iSymOp<vecScrews.size(); ++iSymOp)
		{
			if(bShowMatrices)
				listSymOps->addItem(xtl::print_matrix(vecTrafos[vecScrews[iSymOp]]).c_str());
			else
				listSymOps->addItem(xtl::get_trafo_desc(vecTrafos[vecScrews[iSymOp]]).c_str());
		}
	}
*/

	RecalcBragg();
	CalcTrafo();
}

void SgListDlg::RecalcBragg()
{
	const QListWidgetItem* pItem = listSGs->currentItem();
	if(!pItem) return;

	const int h = spinH->value();
	const int k = spinK->value();
	const int l = spinL->value();
	const bool bOnlyCentring = 0;

	bool (xtl::SpaceGroup<t_real>::*pAllowedFkt)(int, int, int, std::size_t*) const =
		bOnlyCentring ? &xtl::SpaceGroup<t_real>::HasGenReflection
			: &xtl::SpaceGroup<t_real>::HasReflection;

	std::shared_ptr<const xtl::SpaceGroups<t_real>> sgs = xtl::SpaceGroups<t_real>::GetInstance();
	const xtl::SpaceGroups<t_real>::t_vecSpaceGroups* pvecSG = sgs->get_space_groups_vec();
	const unsigned int iSG = pItem->data(Qt::UserRole).toUInt();
	if(iSG >= pvecSG->size())
		return;

	const xtl::SpaceGroup<t_real>* psg = pvecSG->at(iSG);
	std::size_t idxTrafo = 0;
	const bool bForbidden = !(psg->*pAllowedFkt)(h,k,l, &idxTrafo);

	if(bForbidden)
	{
		// mark the symmetry operation responsible for the absence
		for(int iItem=0; iItem<listSGs->count(); ++iItem)
		{
			QListWidgetItem *pItem = listSymOps->item(iItem);
			unsigned int iTrafo = pItem->data(Qt::UserRole).toUInt();
			if(iTrafo == idxTrafo)
			{
				listSymOps->setCurrentItem(pItem);
				break;
			}
		}
	}

	QFont font = spinH->font();
	font.setStrikeOut(bForbidden);
	for(QSpinBox* pSpin : {spinH, spinK, spinL})
		pSpin->setFont(font);
}

void SgListDlg::SearchSG(const QString& qstr)
{
	QList<QListWidgetItem*> lstItems = listSGs->findItems(qstr, Qt::MatchContains);
	if(lstItems.size())
		listSGs->setCurrentItem(lstItems[0], QItemSelectionModel::SelectCurrent);
}


void SgListDlg::CalcTrafo()
{
	listTrafo->clear();

	const QListWidgetItem* pItem = listSGs->currentItem();
	if(!pItem)
		return;

	std::shared_ptr<const xtl::SpaceGroups<t_real>> sgs = xtl::SpaceGroups<t_real>::GetInstance();
	const xtl::SpaceGroups<t_real>::t_vecSpaceGroups* pvecSG = sgs->get_space_groups_vec();
	const unsigned int iSG = pItem->data(Qt::UserRole).toUInt();
	if(iSG >= pvecSG->size())
		return;
	const xtl::SpaceGroup<t_real>* psg = pvecSG->at(iSG);

	xtl::SpaceGroup<t_real>::t_vec vecIn =
		tl::make_vec({spinX->value(), spinY->value(), spinZ->value(), spinW->value()});

	const std::vector<xtl::SpaceGroup<t_real>::t_mat>& vecTrafos = psg->GetTrafos();
	std::vector<xtl::SpaceGroup<t_real>::t_vec> vecUnique;

	listTrafo->addItem(create_header_item("All Transformation Results"));
	for(unsigned int iMat=0; iMat<vecTrafos.size(); ++iMat)
	{
		const xtl::SpaceGroup<t_real>::t_mat& mat = vecTrafos[iMat];
		xtl::SpaceGroup<t_real>::t_vec vec = ublas::prod(mat, vecIn);
		QListWidgetItem* pItem = new QListWidgetItem(xtl::print_vector(vec).c_str());
		pItem->setData(Qt::UserRole, iMat);
		listTrafo->addItem(pItem);

		if(!xtl::is_vec_in_container<std::vector, xtl::SpaceGroup<t_real>::t_vec>(vecUnique, vec))
			vecUnique.push_back(vec);
	}

	listTrafo->addItem(create_header_item("Unique Transformation Results"));
	for(const xtl::SpaceGroup<t_real>::t_vec& vec : vecUnique)
	{
		listTrafo->addItem(xtl::print_vector(vec).c_str());
	}
}


#include "moc_SgListDlg.cpp"
