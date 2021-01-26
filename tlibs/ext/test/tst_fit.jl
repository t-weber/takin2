#
# tlibs Julia module test
# @author Tobias Weber <tobias.weber@tum.de>
# @date 23-apr-2017
# @license GPLv2 or GPLv3
#

include("tlibs.jl")


# load data file
a = tl.loadinstr("/home/tweber/Measurements/mira-mgv2o4-17/data/11797_00025700.dat")
cols = a[1]
data = a[2]
#println(cols)
#println(typeof(data))

idxE = findfirst(cols, "E")
idxCtr = findfirst(cols, "ctr1")

Es = data[:,idxE]
cts = data[:,idxCtr]
cts_err = sqrt.(cts)

println(lpad(cols[idxE], 12), " ", lpad(cols[idxCtr], 12))
for (E, ct) in zip(Es, cts)
	println(lpad(E, 12), " ", lpad(ct, 12))
end
println()


# fit gaussian
fitresult = tl.fit(tl.gauss_model_amp, Es, cts, cts_err, fixed = ["offs"],
	values = Dict("x0" => -1.5, "sigma" => 0.5, "amp" => 100., "offs" => 50.),
	errors = Dict("x0" => 0.25, "sigma" => 0.25, "amp" => 20., "offs" => 10.))
println(fitresult)
println()

data_fine = tl.eval_func(fitresult, minimum(Es), maximum(Es))


# write results
f = open("tst_pts.dat", "w")
for (x,y,yerr) in zip(Es, cts, cts_err)
	println(f, rpad(x, 20), " ", rpad(y, 20), " ", rpad(yerr, 20))
end
close(f)

f = open("tst_fit.dat", "w")
for (x,y) in zip(data_fine[1], data_fine[2])
	println(f, rpad(x, 20), " ", rpad(y, 20))
end
close(f)


# plot results
run(`gnuplot -p -e "plot \"tst_pts.dat\" w yerrorbars pt 7 lc rgb \"#ff0000\",
	\"tst_fit.dat\" with lines lw 2 lc rgb \"#0000ff\""`)
