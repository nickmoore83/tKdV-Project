#= Calculate the transfer function:
For given upstream theta, enforce the matching condition to 
compute the downstream theta.=#

include("math_routines.jl")
include("io_routines.jl")

function transfun(randfile::AbstractString, lamfac::Int)
	rr, nmodes = load(randfile, "rr", "nmodes")
	# Parameters
	eps0 = 0.017
	del0 = 0.23
	D0 = 0.24
	# Constants
	C2,C3 = C2C3(eps0,del0,lamfac)
	cc = ConstantList(C2,C3,D0)
	# Compute the list of Hamiltonians using precomputed H2 and H3
	hdn = Float32.(hamdn(rr,cc))
	hup = Float32.(hamup(rr,cc))

	# Specify list of upstream thetas.
	thup = 0.: 0.5 : 30.
	nth = length(thup)
	thdn = zeros(nth); skup = zeros(nth); skdn = zeros(nth)
	guess = thup[1]
	for nn = 1:2
		meandiff(thdn) = meanham(hdn, hdn, thdn) - meanham(hdn, hup, thup[nn])
		thdn[nn] = find_zero(meandiff, guess, Order1())
		guess = thdn[nn]
		nn>1 ? guess = 2*thdn[nn] - thdn[nn-1] : 0.	# Linear extrapolation
		println("Iteration ", nn, "; thup = ", thup[nn], "; thdn = ", sig(thdn[nn],3))

		# Compute the skewness of eta
		skup[nn] = skewu(rr.H3, hup, thup[nn])
		skdn[nn] = skewu(rr.H3, hdn, thdn[nn])
		println("Upstream skewness = ", sig(skup[nn],3))
		println("Downstream skewness = ", sig(skdn[nn],3), "\n")
	end
	savefile = string("thvars-", string(nmodes), "-", string(lamfac), ".jld")
	save(savefile, "thup", thup, "thdn", thdn, 
		"skup", skup, "skdn", skdn, "nmodes", nmodes, "lamfac", lamfac)
end

function plotstuff(datafile::AbstractString)
	thup,thdn,skup,skdn = load(datafile, "thup", "thdn", "skup", "skdn")
	nmodes,lamfac = load(datafile, "nmodes", "lamfac")
	ratio = skdn./skup
	# Make plots
	mytitle = string("Lambda = ",nmodes, ", N = ", lamfac)
	p1 = plot(thup, thdn, shape = :circle, title = mytitle,
		xlabel = L"\theta^-", ylabel = L"\theta^+", ylims = (0,20))
	p2 = plot(thup, [skup,skdn], title = mytitle,
		xlabel = L"\theta^-", ylabel = "skewness", ylims = (0,2),
		label = ["incoming", "outgoing"], shape = :circle)
	p3 = plot(thup[2:end], ratio[2:end], title = mytitle,
		xlabel = L"\theta^-", ylabel = "skewness ratio", ylims = (0,20),
		label = ["ratio"], shape = :circle)
	savefig(p1,"theta.pdf")
	savefig(p2,"skew.pdf")
	savefig(p3,"ratio.pdf")
	#writedata([thup; thdn], "theta.txt")
end

function output_text(datafile::AbstractString)
	thup,thdn,skup,skdn = load(datafile, "thup", "thdn", "skup", "skdn")
	nmodes,lamfac = load(datafile, "nmodes", "lamfac")
	nthvals = length(thup)
	ratio = skdn./skup
	outdata = [nmodes; lamfac; nthvals; thup; thdn; skup; skdn]
	outfile = string(datafile,".txt")
	writedata(outdata, outfile)
end


## RUN STUFF
#@time( uniform_sample(16, 2, false) )

#transfun("rand-15-2000.jld", 8)

#output_text("thvars-15-8-2000.jld")

