/**
 * popovici calculation
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013-2016
 * @license GPLv2
 *
 * @desc This is a reimplementation in C++ of the file rc_popma.m of the
 *		rescal5 package by Zinkin, McMorrow, Tennant, Farhi, and Wildes (ca. 1995-2007):
 *		http://www.ill.eu/en/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab/
 * @desc see: [pop75] M. Popovici, Acta Cryst. A 31, 507 (1975), doi: 10.1107/S0567739475001088
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

#ifndef __TAKIN_POP_H__
#define __TAKIN_POP_H__

#include "cn.h"

/**
 * TAS parameters in fwhm
 */
struct PopParams : public CNParams
{
	tl::t_length_si<t_real_reso> mono_w;
	tl::t_length_si<t_real_reso> mono_h;
	tl::t_length_si<t_real_reso> mono_thick;
	tl::t_length_si<t_real_reso> mono_curvh;
	tl::t_length_si<t_real_reso> mono_curvv;
	bool bMonoIsCurvedH=0, bMonoIsCurvedV=0;
	bool bMonoIsOptimallyCurvedH=0, bMonoIsOptimallyCurvedV=0;
	unsigned int mono_numtiles_v, mono_numtiles_h;

	tl::t_length_si<t_real_reso> ana_w;
	tl::t_length_si<t_real_reso> ana_h;
	tl::t_length_si<t_real_reso> ana_thick;
	tl::t_length_si<t_real_reso> ana_curvh;
	tl::t_length_si<t_real_reso> ana_curvv;
	bool bAnaIsCurvedH=0, bAnaIsCurvedV=0;
	bool bAnaIsOptimallyCurvedH=0, bAnaIsOptimallyCurvedV=0;
	unsigned int ana_numtiles_v, ana_numtiles_h;

	bool bSampleCub=1;
	tl::t_length_si<t_real_reso> sample_w_q;
	tl::t_length_si<t_real_reso> sample_w_perpq;
	tl::t_length_si<t_real_reso> sample_h;

	bool bSrcRect=1;
	tl::t_length_si<t_real_reso> src_w;
	tl::t_length_si<t_real_reso> src_h;

	bool bDetRect=1;
	tl::t_length_si<t_real_reso> det_w;
	tl::t_length_si<t_real_reso> det_h;

	bool bGuide=0;
	tl::t_angle_si<t_real_reso> guide_div_h;
	tl::t_angle_si<t_real_reso> guide_div_v;

	tl::t_length_si<t_real_reso> dist_mono_sample;
	tl::t_length_si<t_real_reso> dist_sample_ana;
	tl::t_length_si<t_real_reso> dist_ana_det;
	tl::t_length_si<t_real_reso> dist_src_mono;

	tl::t_length_si<t_real_reso> monitor_w;
	tl::t_length_si<t_real_reso> monitor_h;
	tl::t_length_si<t_real_reso> dist_mono_monitor;
};


extern ResoResults calc_pop(const PopParams& pop);
extern ResoResults calc_pop_cn(const CNParams& pop);


// add sample mosaic or other uncertainty to covariance matrix
extern void add_cov_variance(ublas::matrix<t_real_reso>& cov, t_real_reso& R0,
	t_real_reso fwhm_Qx, t_real_reso fwhm_Qy, t_real_reso fwhm_Qz, t_real_reso fwhm_E,
	bool calc_R0 = true);


#endif
