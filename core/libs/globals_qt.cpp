/**
 * globals
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 20-mar-2015
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

#include "globals.h"
#include "globals_qt.h"
#include "tlibs/log/log.h"
#include "tlibs/file/file.h"


QFont g_fontGen("DejaVu Sans",10);
QFont g_fontGfx("DejaVu Sans",10);
//QFont g_fontGL("DejaVu Sans Mono",10);
std::string g_strFontGL("res/fonts/DejaVuSansMono.ttf");
int g_iFontGLSize = 24;


QIcon load_icon(const std::string& strIcon)
{
	std::string strFile = find_resource(strIcon);
	if(strFile != "")
		return QIcon(strFile.c_str());

	return QIcon();
}


QPixmap load_pixmap(const std::string& strIcon)
{
	std::string strFile = find_resource(strIcon);
	if(strFile != "")
		return QPixmap(strFile.c_str());

	return QPixmap();
}
