#
# structure factor test
# @author Tobias Weber <tobias.weber@tum.de>
# @date nov-2017
# @license GPLv2
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

include("sgs.jl")
include("slens.jl")


function structfact(atoms, G)
	F = 0 + 0im

	for atom in atoms
		F += exp(2.*pi * 1im * dot(G, atom))
	end

	return F
end


function has_pos(allpos, pos)
	for curpos in allpos
		if(curpos == pos)
			return true
		end
	end

	return false
end


function generate_pos(pos, trafos)
	allpos = []

	for trafo in trafos
		rot = trafo[1:3, 1:3]
		trans = trafo[1:3, 4:4]

		newpos = rot * pos + trans

#		for idx in range(1, length(newpos))
#			while(newpos[idx] < 0.)
#				newpos[idx] += 1.
#			end
#			while(newpos[idx] > 1.)
#				newpos[idx] -= 1.
#			end
#		end

		if !has_pos(allpos, newpos)
			push!(allpos, newpos)
		end
	end

	#println(allpos)
	return allpos
end


atom_1 = "Mn"
atom_2 = "Si"
atompos_1 = [ 0.14; 0.14; 0.14 ]
atompos_2 = [ 0.85; 0.85; 0.85 ]
sgroup = "P2_13"

slen_1 = slens[atom_1]["coh"]
slen_2 = slens[atom_2]["coh"]
trafos = sgs[sgroup]["trafos"]

allatompos_1 = generate_pos(atompos_1, trafos)
allatompos_2 = generate_pos(atompos_2, trafos)


G = [ 1.; 1.; 0. ]

F_1 = structfact(allatompos_1, G) * slen_1
F_2 = structfact(allatompos_2, G) * slen_2
F = F_1 + F_2


println(atom_1, " atoms: ", allatompos_1)
println(atom_2, " atoms: ", allatompos_2)

println("G = ", G, ", F = ", F, ", |F| = ", sqrt(F*conj(F)))
