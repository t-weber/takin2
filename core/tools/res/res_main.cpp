/**
 * reso tool
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013, 2014
 * @copyright GPLv2
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
#include "ResoDlg.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/log/log.h"

#ifdef Q_WS_X11
	extern "C" int XInitThreads();
#endif

int main(int argc, char** argv)
{
	try
	{
		tl::log_info("Starting up resolution tool.");

		#ifdef Q_WS_X11
			XInitThreads();
		#endif

		tl::init_spec_chars();
		QApplication app(argc, argv);

		std::setlocale(LC_ALL, "C");
		QLocale::setDefault(QLocale::English);

		ResoDlg dlg(0);
		dlg.show();
		int iRet = app.exec();

		tl::deinit_spec_chars();
		tl::log_info("Shutting down resolution tool.");
		return iRet;
	}
	catch(const std::exception& ex)
	{
		tl::log_crit(ex.what());
	}

	return -1;
}
