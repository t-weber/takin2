#
# tlibs Julia module test
# @author Tobias Weber <tobias.weber@tum.de>
# @date 23-apr-2017
# @license GPLv2 or GPLv3
#
# ----------------------------------------------------------------------------
# tlibs -- a physical-mathematical C++ template library
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) version 3.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------------
#

include("tlibs.jl")


# a*x^2 + b*x + c
# 2*a*x + b = 0  =>  x_extr = -b/(2a)
function fkt(x, y)
	return 0.5*x^2. + 10.*x + y
end


minresult = tl.minimise(fkt, fixed = ["y"],
	values = Dict("x" => 0., "y" => 15.),
	errors = Dict("x" => 10., "y" => 0.))

println(minresult)

