/**
 * implementation of the eckold-sobolev algo
 *
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2015
 * @license GPLv2
 *
 * @desc for algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014), doi: 10.1016/j.nima.2014.03.019
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

#ifndef __TAKIN_ECK_H__
#define __TAKIN_ECK_H__

#include "pop.h"

/**
 * TAS parameters in fwhm
 */
struct EckParams : public PopParams
{
	tl::t_length_si<t_real_reso> pos_x, pos_y, pos_z;

	// vertical scattering in k_f
	bool bKfVertical = false;
};


extern ResoResults calc_eck(const EckParams& eck);


#endif
