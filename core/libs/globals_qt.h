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

#ifndef __TAKIN_GLOBALS_QT_H__
#define __TAKIN_GLOBALS_QT_H__

#include <string>

#include <QIcon>
#include <QPixmap>
#include <QFont>


extern QIcon load_icon(const std::string& strIcon);
extern QPixmap load_pixmap(const std::string& strPixmap);

extern QFont g_fontGen, g_fontGfx/*, g_fontGL*/;
extern std::string g_strFontGL;
extern int g_iFontGLSize;


#endif
