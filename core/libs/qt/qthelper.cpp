/**
 * qt helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2016
 * @license GPLv2
 */

#include "qthelper.h"
//#include "globals.h"
//#include "globals_qt.h"
#include "tlibs/math/math.h"
#include "tlibs/string/string.h"
#include "tlibs/helper/misc.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <memory>

#include <QFileDialog>
#include <QMessageBox>
#include <QMouseEvent>


bool save_table(const char* pcFile, const QTableWidget* pTable)
{
	std::ofstream ofstr(pcFile);
	if(!ofstr)
		return false;

	const int iNumCols = pTable->columnCount();
	const int iNumRows = pTable->rowCount();

	// item lengths
	std::unique_ptr<int[]> ptrMaxTxtLen(new int[iNumCols]);
	for(int iCol=0; iCol<iNumCols; ++iCol)
		ptrMaxTxtLen[iCol] = 0;

	for(int iCol=0; iCol<iNumCols; ++iCol)
	{
		const QTableWidgetItem *pItem = pTable->horizontalHeaderItem(iCol);
		ptrMaxTxtLen[iCol] = std::max(pItem ? pItem->text().length() : 0, ptrMaxTxtLen[iCol]);
	}

	for(int iRow=0; iRow<iNumRows; ++iRow)
	{
		for(int iCol=0; iCol<iNumCols; ++iCol)
		{
			const QTableWidgetItem *pItem = pTable->item(iRow, iCol);
			ptrMaxTxtLen[iCol] = std::max(pItem ? pItem->text().length() : 0, ptrMaxTxtLen[iCol]);
		}
	}


	// write items
	for(int iCol=0; iCol<iNumCols; ++iCol)
	{
		const QTableWidgetItem *pItem = pTable->horizontalHeaderItem(iCol);
		ofstr << std::setw(ptrMaxTxtLen[iCol]+4) << (pItem ? pItem->text().toStdString() : "");
	}
	ofstr << "\n";

	for(int iRow=0; iRow<iNumRows; ++iRow)
	{
		for(int iCol=0; iCol<iNumCols; ++iCol)
		{
			const QTableWidgetItem *pItem = pTable->item(iRow, iCol);
			ofstr << std::setw(ptrMaxTxtLen[iCol]+4) << (pItem ? pItem->text().toStdString() : "");
		}

		ofstr << "\n";
	}

	return true;
}


// ----------------------------------------------------------------------------


#include <QStandardPaths>

std::vector<std::string> get_qt_std_path(QtStdPath path)
{
	QStandardPaths::StandardLocation iLoc;
	switch(path)
	{
		case QtStdPath::FONTS:
			iLoc = QStandardPaths::FontsLocation;
			break;
		case QtStdPath::HOME:
		default:
			iLoc = QStandardPaths::HomeLocation;
			break;
	}

	QStringList lst = QStandardPaths::standardLocations(iLoc);

	std::vector<std::string> vecPaths;
	for(int iStr=0; iStr<lst.size(); ++iStr)
		vecPaths.push_back(lst.at(iStr).toStdString());
	return vecPaths;
}


// ----------------------------------------------------------------------------


void focus_dlg(QDialog* pDlg)
{
	if(!pDlg) return;

	pDlg->show();
	pDlg->raise();
	pDlg->activateWindow();
}
