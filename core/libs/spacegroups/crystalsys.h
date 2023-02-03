/**
 * Crystal Systems
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

#ifndef __XTL_CRYSSYS_H__
#define __XTL_CRYSSYS_H__

namespace xtl {

enum CrystalSystem
{
	CRYS_NOT_SET,

	CRYS_TRICLINIC,		// all params free
	CRYS_MONOCLINIC,	// beta=gamma=90
	CRYS_ORTHORHOMBIC,	// alpha=beta=gamma=90
	CRYS_TETRAGONAL,	// a=b, alpha=beta=gamma=90
	CRYS_TRIGONAL,		// a=b=c, alpha=beta=gamma
	CRYS_HEXAGONAL,		// a=b, gamma=120, alpha=beta=90
	CRYS_CUBIC		// a=b=c, alpha=beta=gamma=90
};


extern const unsigned int* get_crystal_system_start_indices();
extern const char** get_crystal_system_names(bool bCapital=0);
extern CrystalSystem get_crystal_system_from_space_group(unsigned int iSGNr);
extern CrystalSystem get_crystal_system_from_laue_group(const char* pcLaue);
extern const char* get_crystal_system_name(CrystalSystem ty);

}

#ifndef NO_QT
#include <QLineEdit>

namespace xtl {

extern void set_crystal_system_edits(CrystalSystem crystalsys,
	QLineEdit* pA, QLineEdit* pB, QLineEdit* pC,
	QLineEdit* pAlpha, QLineEdit* pBeta, QLineEdit* pGamma,
	QLineEdit* pAReci=nullptr, QLineEdit* pBReci=nullptr, QLineEdit* pCReci=nullptr,
	QLineEdit* pAlphaReci=nullptr, QLineEdit* pBetaReci=nullptr, QLineEdit* pGammaReci=nullptr);
}
#endif

#endif
