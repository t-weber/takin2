/**
 * Log Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 14-oct-2017
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

#include <QScrollBar>
#include "LogDlg.h"
#include "tlibs/file/file.h"
#include <fstream>


LogDlg::LogDlg(QWidget* pParent, QSettings *pSett, const std::string& strLogFile)
	: QDialog(pParent), m_pSettings(pSett)
{
	this->setupUi(this);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);

		if(m_pSettings->contains("log/geo"))
			restoreGeometry(m_pSettings->value("log/geo").toByteArray());
	}

	if(tl::file_exists(strLogFile.c_str()))
	{
		m_pFileWatcher.reset(new QFileSystemWatcher(this));
		m_pFileWatcher->addPath(strLogFile.c_str());

		QObject::connect(m_pFileWatcher.get(), &QFileSystemWatcher::fileChanged,
			this, &LogDlg::LogFileChanged);

		LogFileChanged(strLogFile.c_str());
	}
	else
	{
		editOut->setPlainText("Error: No log data was found!");
	}
}


/**
 * log file changed
 */
void LogDlg::LogFileChanged(const QString& _strLogFile)
{
	std::string strLogFile = _strLogFile.toStdString();
	std::ifstream ifstr(strLogFile);
	if(!ifstr)
		return;

	std::size_t iSize = std::size_t(tl::get_file_size(ifstr));
	std::unique_ptr<char[]> pcLog(new char[iSize+1]);

	ifstr.read((char*)pcLog.get(), iSize);
	pcLog[iSize] = 0;
	editOut->setPlainText(pcLog.get());

	// scroll to end
	QScrollBar *pSB = editOut->verticalScrollBar();
	pSB->setValue(pSB->maximum());
}


void LogDlg::accept()
{
	if(m_pSettings)
		m_pSettings->setValue("log/geo", saveGeometry());

	QDialog::accept();
}

#include "moc_LogDlg.cpp"
