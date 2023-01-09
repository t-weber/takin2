/**
 * Scan viewer
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date mar-2015
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

#include <clocale>
#include <QLocale>
#include <QApplication>

#include "scanviewer.h"
#include "tlibs/math/rand.h"
#include "libs/version.h"


int main(int argc, char** argv)
{
	tl::init_rand();

	QApplication app(argc, argv);

	app.setApplicationName("Takin/Scanviewer");
	app.setApplicationVersion(TAKIN_VER);

	std::setlocale(LC_ALL, "C");
	QLocale::setDefault(QLocale::English);

	ScanViewerDlg dlg(0);
	dlg.setWindowFlags(Qt::Window);
	dlg.show();

	int iRet = app.exec();
	return iRet;
}
