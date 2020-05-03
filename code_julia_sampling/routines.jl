#= All routines assume that the signal uu is real and has zero mean.
Therefore the Fourier transform, uhat, only requires k=1,...,nmodes. 
The negative modes are complex conjugates and the k=0 mode vanishes.=#

#---------- Basic Stuff ----------#
using Distributions
using Roots
using DelimitedFiles
using LinearAlgebra
using FFTW
using JLD
using Plots
using LaTeXStrings
# Set the data folder.
data_folder() = "./Data"
# Round to a set number of significant digits
sig(var, sigdig::Int) = round(var,sigdigits=sigdig)
#---------------------------------#

#---------- Real FFT Routines ----------#
# Note: Optimally, the length of uhat and uu should be a power of 2.
# Note on symint: to use the symmetric interval [-pi,pi], odd modes must be negated.
# Transform from physical to spectral space; only used in benchmark of H3.
function realfft(uu::Vector{Float64}, symint::Bool=true)
	uhat = rfft(uu)/length(uu)
	@assert(abs(uhat[1])/maximum(abs,uhat) < 1e-6)
	uhat[end] *= 0.5
	uhat = uhat[2:end]
	if(symint) uhat[1:2:end] .*= -1 end
	return uhat
end
# Transform from spectral to physical space; used in ham3fft.
# Note: uh must be created to prevent input uhat from being destroyed.
function irealfft(uhat::Vector{ComplexF64}, symint::Bool=true)
	uh = [0; uhat[1:end-1]; 2*uhat[end]]
	if(symint) uh[2:2:end] .*= -1 end
	uu = irfft(uh, 2*length(uhat))	# Note: irfft destroys input uh.
	return uu*length(uu)
end
# Upsampled version of irealfft.
function ifftup(uhat::Vector{ComplexF64}, symint::Bool=true) 
	return irealfft([uhat; zeros(length(uhat))], symint)
end
#---------------------------------------#

#---------- Basic Hamiltonian Routines ----------#
# Compute the energy, E = 1/2 int u^2 dx.
energy(uhat::Vector{ComplexF64}) = 2*pi*norm(uhat)^2
# Compute H2 = 1/2 int u_x^2 dx
ham2(uhat::Vector{ComplexF64}) = energy(im*(1:length(uhat)).*uhat)
# Compute H3 using FFT to physical space, H3 = 1/6 int u^3 dx.
function ham3fft(uhat::Vector{ComplexF64})
	uu = ifftup(uhat)
	return 1/6 * sum(uu.^3) * 2*pi/length(uu)
end
# Recursive computation of H3 in spectral space.
function ham3rec(uhat::Vector{ComplexF64})
	# If only one mode, then H3 is zero.
	length(uhat) == 1 && return 0.0
	# Otherwise use recursion.
	H3 = ham3rec(uhat[1:end-1])
	usum = sum( uhat[nn]*uhat[end-nn] for nn=1:lastindex(uhat)-1 )
	H3 += 2*pi*real( conj(uhat[end]) *  usum)
	return H3
end
# Choose which algorithm to use for H3.
ham3(uhat::Vector{ComplexF64}) = ham3rec(uhat) 
#---------------------------------------#

#---------- Sampling Routines ----------#
# Sampled state
struct RandList; rvar::Array; H2::Vector; H3::Vector; end
# Get uhat from the array of random values.
function getuhat(rvar::Array{Float64}, nn::Int)
	uhat = rvar[:,1,nn] + im*rvar[:,2,nn]
	return uhat/sqrt(energy(uhat))
end
# Sample a single sweep from a uniform distribution on the hypershpere E=1.
function sample_one_sweep(nmodes::Int, nsamp::Int, zerolast::Bool)
	rvar = randn(nmodes,2,nsamp)
	# Zero out the last mode as is done in Matlab DNS.
	if zerolast
		rvar[nmodes,:,:] *= 0.
	end
	# Compute H3 and H2 for each sample in a parallel for loop.
	H2vec, H3vec = [zeros(Float64,nsamp) for nn=1:2]
	## Later parralelize this step. #@parallel
	for nn = 1:nsamp
		uhat = getuhat(rvar,nn)
		H2vec[nn] = ham2(uhat)
		H3vec[nn] = ham3(uhat)
	end
	rlist = RandList(Float32.(rvar), H2vec, H3vec)
	return rlist
end
# Sample several sweeps for better memory usage.
function sample_many_sweeps(nmodes::Int, nsweeps::Int, zerolast::Bool)
	samp_per = 10^5
	totsamp = nsweeps*samp_per
	H2vec, H3vec = [zeros(Float32,totsamp) for nn=1:2]
	for nn = 1:nsweeps
		rlist = sample_one_sweep(nmodes, samp_per, zerolast)
		idx_range = (nn-1)*samp_per+1 : nn*samp_per
		H2vec[idx_range] = Float32.(rlist.H2)
		H3vec[idx_range] = Float32.(rlist.H3)
		println("Sweep ",nn,", ",sig(100*nn/nsweeps,3),"% completed.")
	end
	return H2vec, H3vec, totsamp
end
#= Sample from a uniform distribution by making several sweeps.
Note: zerolast zeroes the last mode as is done in the Matlab DNS. =#
function uniform_sample(nmodes::Int, nsweeps::Int; zerolast::Bool=true)
	(H2vec, H3vec, totsamp), cputime = @timed( sample_many_sweeps(nmodes, nsweeps, zerolast) )
	println("\nCompleted sampling with nmodes = ", nmodes, ", nsweeps = ", nsweeps)
	println("CPU time = ", sig(cputime/60,3), " minutes.\n\n")
	zerolast ? zstr = "z" : zstr = ""
	savefile = string(data_folder(),"rand-",string(nmodes),zstr,"-",string(nsweeps),".jld")
	save(savefile, "rr", RandList([],H2vec,H3vec), 
		"nmodes", nmodes, "nsweeps", nsweeps, "totsamp", totsamp, "cputime", cputime)
end
# OLD RELEVANT CODE SNIPPET FOR SAVEMICRO:
#savemicro ? rsave = Float32.(rvar) : rsave = []
#---------------------------------------#


#---------- Highl-level Hamiltonian Routines ----------#
# List of constants
struct ConstList; C2::Float64; C3::Float64; D0::Float64; end
#= Compute the constanst C2 and C3. =#
function C2C3(eps0::Float64, del0::Float64, lamfac::Int)
	C2 = (2/3) * pi^2 * del0/lamfac^2
	C3 = 1.5 * sqrt(pi) * eps0/del0
	return C2,C3
end
# Compute the Hamiltonian. H3 and H2 are real and can be numbers or vectors.
hamiltonian(H2, H3, C2, C3, D0) = C2*D0^(1/2)*H2 - C3*D0^(-3/2)*H3
# Compute the upstream and downstream Hamiltonians.
hamup(rr::RandList, cc::ConstList) = hamiltonian(rr.H2, rr.H3, cc.C2, cc.C3, 1.)
hamdn(rr::RandList, cc::ConstList) = hamiltonian(rr.H2, rr.H3, cc.C2, cc.C3, cc.D0)
# Mean of ham under Gibbs measure hamgibbs.
function gibbs_mean(quantity, hamgibbs, beta::Float64)
	return dot(exp.(-beta*hamgibbs), quantity) / sum(exp.(-beta*hamgibbs))
end
# Compute the skewness of the displacement, u or eta.
function skewu(H3, hamgibbs, beta::Float64)
	return 3*pi^(1/2)*gibbs_mean(H3, hamgibbs, beta)
end
#---------------------------------------#

#---------- Highest-level Routines ----------#
#= Calculate the transfer function:
Given an upstream beta, enforce the matching condition to compute the downstream beta.=#
function transfun(randfile::AbstractString, lamfac::Int)
	# Set Parameters.
	thup = 0.: 0.5 : 30.
	eps0 = 0.017
	del0 = 0.23
	D0 = 0.24
	# Load data.
	randfile = string(data_folder(), randfile)
	rr, nmodes, nsweeps = load(randfile, "rr", "nmodes", "nsweeps")
	# Set Consts.
	C2,C3 = C2C3(eps0,del0,lamfac)
	cc = ConstList(C2,C3,D0)
	# Compute the list of Hamiltonians using precomputed H2 and H3.
	hdn = Float32.(hamdn(rr,cc))
	hup = Float32.(hamup(rr,cc))
	# Initialize variables.
	nth = length(thup)
	thdn, skup, skdn = [zeros(nth) for idx=1:3]
	guess = thup[1]
	# Loop over the beta values to compute the transfer function.
	cputime = @elapsed(
	for nn = 1:nth
		meandiff(thdn) = gibbs_mean(hdn, hdn, thdn) - gibbs_mean(hdn, hup, thup[nn])
		thdn[nn] = find_zero(meandiff, guess, Order1())
		println("Step ", nn, ": thup = ", thup[nn], ", thdn = ", sig(thdn[nn],3))
		# Use linear extrapolation to make the next guess.
		nn>1 ? guess = 2*thdn[nn] - thdn[nn-1] : guess = thdn[nn]
		# Compute the skewness of eta
		skup[nn] = skewu(rr.H3, hup, thup[nn])
		skdn[nn] = skewu(rr.H3, hdn, thdn[nn])
		println("Upstream skewness = ", sig(skup[nn],3))
		println("Downstream skewness = ", sig(skdn[nn],3), "\n")
	end)
	println("\nCompleted transfer with nmodes = ", nmodes, ", nsweeps = ", nsweeps, " lamfac = ", lamfac)
	println("CPU time = ", sig(cputime/60,3), " minutes.\n\n")	
	savefile = string(data_folder(), "thvars-", string(nmodes), "z-", 
		string(lamfac), "-", string(nsweeps), ".jld")
	save(savefile, "thup", thup, "thdn", thdn, "skup", skup, "skdn", skdn, 
		"nmodes", nmodes, "lamfac", lamfac, "nsweeps", nsweeps, "cputime", cputime)
end

# Output the data to a text file.
function output_text(datafile::AbstractString)
	datafile = string(data_folder(), datafile)
	thup,thdn,skup,skdn = load(datafile, "thup", "thdn", "skup", "skdn")
	nmodes,lamfac,nsweeps = load(datafile, "nmodes", "lamfac", "nsweeps")
	nthvals = length(thup)
	ratio = skdn./skup
	outdata = [nmodes; lamfac; nthvals; thup; thdn; skup; skdn]
	outfile = string(data_folder(), "transf-", string(nmodes), "z-", 
		string(lamfac), "-", string(nsweeps), ".txt")
	writedlm(outfile, outdata)
end
