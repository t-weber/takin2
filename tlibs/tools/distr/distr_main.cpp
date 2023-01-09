/**
 * random distributions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2016
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#include "distr.h"
#include "log/log.h"
#include <memory>

int main(int argc, char** argv)
{
	try
	{
		std::unique_ptr<QApplication> pApp = std::make_unique<QApplication>(argc, argv);
		::setlocale(LC_ALL, "C");
		QLocale::setDefault(QLocale::English);

		std::unique_ptr<DistrDlg> pDlg = std::make_unique<DistrDlg>();
		pDlg->show();
		return pApp->exec();
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
	return -1;
}
