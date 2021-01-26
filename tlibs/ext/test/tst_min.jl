#
# tlibs Julia module test
# @author Tobias Weber <tobias.weber@tum.de>
# @date 23-apr-2017
# @license GPLv2 or GPLv3
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

