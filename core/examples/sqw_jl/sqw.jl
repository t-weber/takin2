#
# Sample Julia S(q,w) module for ferromagnetic dispersions
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
# @date dec-2016
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

import LinearAlgebra


# from scipy.constants
kB = 0.08617


# example dispersion
function disp_ferro(q, D, offs)
	return D*q^2.0 + offs
end

# Gaussian peak
function gauss(x, x0, sig, amp)
	norm = (sqrt(2.0*pi) * sig)
	return amp * exp(-0.5*((x-x0)/sig)^2.0) / norm
end


# Bose factor
function bose(E, T)
	n = 1.0/(exp(abs(E)/(kB*T)) - 1.0)
	if E >= 0.0
		n += 1.0
	end
	return n
end

# Bose factor which is cut off below Ecut
function bose_cutoff(E, T, Ecut=0.02)
	Ecut = abs(Ecut)

	b = 0.0
	if abs(E) < Ecut
		b = bose(sign(E)*Ecut, T)
	else
		b = bose(E, T)
	end

	return b
end



# -----------------------------------------------------------------------------

#
# global variables which can be accessed / changed by Takin
#
g_G = vec([1.0, 1.0, 0.0])	# Bragg peak
g_D = 50.0			# magnon stiffness
g_offs = 0.0			# energy gap
g_sig = 0.02			# linewidth
g_S0 = 10.0			# intensity

g_inc_sig = 0.02	# incoherent width
g_inc_amp = 10.0	# incoherent intensity

g_T = 100.0		# temperature
g_bose_cut = 0.02	# Bose cutoff



# -----------------------------------------------------------------------------

#
# the init function is called after Takin has changed a global variable (optional)
#
function TakinInit()
	println("Calling TakinInit")
end


#
# dispersion E(Q) and weight factor (optional)
#
function TakinDisp(h::Float64, k::Float64, l::Float64)
	# momentum
	Q = vec([h,k,l])
	# reduced momentum
	q = LinearAlgebra.norm(Q - g_G)

	# energy
	E_peak = disp_ferro(q, g_D, g_offs)
	# weight
	w_peak = 1.0
	return [[E_peak, -E_peak], [w_peak, w_peak]]
end


#
# called for every Monte-Carlo point
#
function TakinSqw(h::Float64, k::Float64, l::Float64, E::Float64)::Float64
	#println("Calling TakinSqw(", h, ", ", k, ", ", l, ", ", E, ") -> ", S)
	Es, ws = TakinDisp(h,k,l)

	S_p = gauss(E, Es[1], g_sig, g_S0*ws[1])
	S_m = gauss(E, Es[2], g_sig, g_S0*ws[2])
	incoh = gauss(E, 0.0, g_inc_sig, g_inc_amp)

	b = 1.0
	#b = bose_cutoff(E, g_T, g_bose_cut)
	S = (S_p + S_m)*b + incoh
	return Float64(S)
end


#
# background function, called for every nominal (Q, E) point (optional)
#
function TakinBackground(h::Float64, k::Float64, l::Float64, E::Float64)::Float64
	#println("Calling TakinBackground(", h, ", ", k, ", ", l, ", ", E, ") -> ", S)
	return Float64(0.)
end



# -----------------------------------------------------------------------------
# test
#
#TakinInit()
#println(TakinDisp(1.0, 1.0, 0.0))
#println(TakinSqw(1.0, 1.0, 0.0, 0.0))
