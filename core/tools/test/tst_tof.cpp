/**
 * @author Tobias Weber <tobias.weber@tum.de>
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

// gcc -DNO_QT -I. -I../.. -o tst_tof ../../tools/test/tst_tof.cpp ../../tools/res/vio.cpp ../../tlibs/log/log.cpp -lstdc++ -std=c++11 -lstdc++ -lm

#include "tools/res/vio.h"
#include <iostream>

int main()
{
	VioParams parms;

	parms.ki = 1.4 / tl::get_one_angstrom<double>();
	parms.kf = 1.4 / tl::get_one_angstrom<double>();
	//parms.E = 0. * tl::get_one_meV<double>();

	parms.len_pulse_mono = 10. * tl::get_one_meter<double>();
	parms.len_mono_sample = 1. * tl::get_one_meter<double>();
	parms.len_sample_det = 5. * tl::get_one_meter<double>();

	parms.twotheta = tl::d2r(75.) * tl::get_one_radian<double>();
	parms.twotheta_i = 0. * tl::get_one_radian<double>();
	parms.angle_outplane_i = 0. * tl::get_one_radian<double>();
	parms.angle_outplane_f = 0. * tl::get_one_radian<double>();

	parms.sig_len_pulse_mono = 0.01 * tl::get_one_meter<double>();
	parms.sig_len_mono_sample = 0.01 * tl::get_one_meter<double>();
	parms.sig_len_sample_det = 0.01 * tl::get_one_meter<double>();

	parms.sig_pulse = 50e-6 * tl::get_one_second<double>();
	parms.sig_mono = 5e-6 * tl::get_one_second<double>();
	parms.sig_det = 5e-6 * tl::get_one_second<double>();

	parms.sig_twotheta_i = tl::m2r(30.) * tl::get_one_radian<double>();
	parms.sig_twotheta_f = tl::m2r(30.) * tl::get_one_radian<double>();
	parms.sig_outplane_i = tl::m2r(30.) * tl::get_one_radian<double>();
	parms.sig_outplane_f = tl::m2r(30.) * tl::get_one_radian<double>();

	ResoResults res = calc_vio(parms);

	return 0;
}
