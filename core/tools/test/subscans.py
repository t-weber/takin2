#
# subtract datasets
#
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
# @date feb-2022
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------
#

import numpy as np
import numpy.linalg as la

import scipy as sp
import scipy.constants as const


# kB in meV/K
kB = const.k / const.e * 1e3


# Bose factor
def bose(E, T):
	n = 1./(np.exp(abs(E)/(kB*T)) - 1.)
	n = np.where(n >= 0., n+1., n)
	return n


# Bose factor which is cut off below Ecut
def bose_cutoff(E, T, Ecut=0.02):
	Ecut = abs(Ecut)

	if abs(E) < Ecut:
		b = bose(np.sign(E)*Ecut, T)
	else:
		b = bose(E, T)

	return b


E_col     = 0
I_col     = 1
Ierr_col  = 2

file_in1  = "/Users/t_weber/Projects/cso/convo/in8/data/38552_m35_m35_05_B110_norm.dat"
#file_in1 = "/Users/t_weber/Projects/cso/convo/in8/data/38539_m35_m35_05_B001_norm.dat"
file_in2  = "/Users/t_weber/Projects/cso/convo/in8/data/38576_m35_m35_05_T80_norm.dat"
#file_in2 = "/Users/t_weber/Projects/cso/convo/in8/data/38582_m35_m35_05_T80_B0_norm.dat"
file_out  = "diff.dat"

T_in1     = 10.
T_in2     = 80.

offs      = 5.5e-4


if __name__ == "__main__":
	dat1 = np.loadtxt(file_in1)
	dat2 = np.loadtxt(file_in2)

	E1 = dat1[:, E_col]
	bose1 = bose(E1, T_in1)
	I1 = dat1[:, I_col] / bose1
	I1_err = dat1[:, Ierr_col] / bose1

	E2 = dat2[:, E_col]
	bose2 = bose(E2, T_in2)
	I2 = dat2[:, I_col] / bose2
	I2_err = dat2[:, Ierr_col] / bose2

	if I1.size != I2.size:
		print("Error: Unequal array sizes.")
		exit(-1)

	E_diff = np.abs(E1 - E2)
	if np.any(E_diff > 0.025):
		print("Error: Cannot subtract non-matching energies.")
		exit(-1)

	I_sub = I1 - I2 + offs
	I_sub_err = np.sqrt(I1_err**2. + I2_err**2.)
	#I_sub = np.where(I_sub >= 0., I_sub, 0.)

	#import matplotlib.pyplot as plt
	#plt.errorbar(E1, I_sub, yerr=I_sub_err, fmt='.')
	#plt.show()

	np.savetxt(file_out, np.array([E1, I_sub, I_sub_err]).T)
