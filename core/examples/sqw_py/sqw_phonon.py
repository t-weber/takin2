#
# Sample Python S(q,w) module for simple acoustic phonons
#
# @author Tobias Weber <tweber@ill.fr>
# @license GPLv2
# @date nov-2019
#

import math as m

import numpy as np
import numpy.linalg as la
from numpy import array	# in global namespace so that Takin can access it

import scipy as sp
import scipy.constants as const



# -----------------------------------------------------------------------------
# dispersion
# -----------------------------------------------------------------------------

# kB in meV/K
kB = const.k / const.e * 1e3


# dispersion relations
def disp_phonon(q, amp, freq, offs):
	return np.abs(amp*np.sin(freq*q)) + offs



# Bose factor
def bose(E, T):
	n = 1./(m.exp(abs(E)/(kB*T)) - 1.)
	if E >= 0.:
		n += 1.
	return n

# Bose factor which is cut off below Ecut
def bose_cutoff(E, T, Ecut=0.02):
	Ecut = abs(Ecut)

	if abs(E) < Ecut:
		b = bose(np.sign(E)*Ecut, T)
	else:
		b = bose(E, T)

	return b



# Gaussian peak
def gauss(x, x0, sig, amp):
	norm = (np.sqrt(2.*m.pi) * sig)
	return amp * np.exp(-0.5*((x-x0)/sig)**2.) / norm


#
# peak shape of a damped harmonic oscillator
# see: B. Fak, B. Dorner, Physica B 234-236 (1997) pp. 1107-1108
#
def DHO(E, T, E0, hwhm, amp):
	return np.abs(bose(E, T)*amp/(E0*np.pi) *
		(hwhm/((E-E0)**2. + hwhm**2.) - hwhm/((E+E0)**2. + hwhm**2.)))

# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Takin interface
# -----------------------------------------------------------------------------

# global variables which can be accessed / changed by Takin
g_G = np.array([4., 4., 0.])	# Bragg peak

g_amp = 20.		# amplitude of sinusoidal dispersion
g_freq = np.pi/2.	# frequency of sinusoidal dispersion
g_offs = 0.		# energy gap
g_HWHM = 0.02		# linewidth
g_S0 = 1.		# intensity

g_inc_sig = 0.02	# incoherent width
g_inc_amp = 1.		# incoherent intensity

g_T = 300.		# temperature

g_bose_cut = 0.02	# cutoff energy for Bose factor


#
# the init function is called after Takin has changed a global variable (optional)
#
def TakinInit():
	print("Init: G=" + repr(g_G) + ", T=" + repr(g_T))


#
# dispersion E(Q) and weight factor (optional)
#
def TakinDisp(h, k, l):
	E_peak = 0.		# energy
	w_peak = 1.		# weight

	try:
		Q = np.array([h,k,l])
		q = la.norm(Q - g_G)
		E_peak = disp_phonon(q, g_amp, g_freq, g_offs)
	except ZeroDivisionError:
		return [0., 0.]

	return [[E_peak, -E_peak], [w_peak, w_peak]]


#
# S(Q,E) function, called for every Monte-Carlo point
#
def TakinSqw(h, k, l, E):
	try:
#		print("h={0}, k={1}, l={2}, E={3}".format(h,k,l,E))

		[Ep_peak, Em_peak], [wp_peak, wm_peak] = TakinDisp(h,k,l)

#		S_p = DHO(E, g_T, Ep_peak, g_HWHM, g_S0*wp_peak)
#		S_m = DHO(E, g_T, Em_peak, g_HWHM, g_S0*wm_peak)
#
#		S = (S_p + S_m)*bose_cutoff(E, g_T, g_bose_cut) + incoh

		S = np.abs(DHO(E, g_T, Ep_peak, g_HWHM, g_S0*wp_peak))
		incoh = gauss(E, 0., g_inc_sig, g_inc_amp)

		S = S + incoh
#		print("S={0}".format(S))
		return S
	except ZeroDivisionError:
		return 0.

# -----------------------------------------------------------------------------


import os
print("Script working directory: " + os.getcwd())

# test
#print(TakinSqw(4.1, 3.9, 0., 0.4))
