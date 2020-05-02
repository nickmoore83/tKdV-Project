#---------- Basic Stuff ----------#
using Distributions
using Roots
using DelimitedFiles
using LinearAlgebra
using FFTW
using JLD
using Plots
using LaTeXStrings
#---------------------------------#

# Set the data folder.
function data_folder()
	return "./Data/"
end
# Round to a set number of significant digits
function sig(var, sigdig)
	return round(var,sigdigits=sigdig)
end

#---------- Real FFT Routines ----------#
#= Compute the FFT between physical and spectral space assuming
the signal u is real and has zero mean. =#
#= Basic realfft to go from physical to spectral space. 
This routine is only used in the benchmark. =#
function realfft(uu::Vector{Float64})
	uhat = rfft(uu)/length(uu)
	@assert(abs(uhat[1])/maximum(abs,uhat) < 1e-6) 
	return uhat[2:end]
end
#= Basic irealfft to go from spectral to physical space. =#
function irealfft(uhat::Vector{Complex{Float64}})
	uu = irfft([0; uhat[:]], 2*length(uhat))
	return uu*length(uu)
end
#= Upsampled version of irealfft. =#
function ifftup(uhat::Vector{Complex{Float64}})
	uhat = [uhat; zeros(eltype(uhat), length(uhat))]
	return irealfft(uhat)
end
#---------------------------------------#

#---------- Hamiltonian Routines ----------#
#= In all routines, the physical signal is assumed to be real with zero momentum. 
Therefore, the Fourier transform only requires k=1,...,nmodes. 
The negative modes are given by the complex conjugates (CC), and the k=0 mode vanishes. =#
#---------- Low-level Hamiltonian Routines ----------#
#= Compute the energy, E = 1/2 int u^2 dx. =#
function energy(uhat::Vector{Complex{Float64}})
	return 2*pi*norm(uhat)^2
end
# Compute H2 = 1/2 int u_x^2 dx
function ham2(uhat::Vector{Complex{Float64}})
	kvec = 1:length(uhat)
	uxhat = im*kvec.*uhat
	return energy(uxhat)
end
#= Compute H3 using FFT to physical space, H3 = 1/6 int u^3 dx. =#
function ham3(uhat::Vector{Complex{Float64}})
	uu = ifftup(uhat)
	return 1/6 * sum(uu.^3) * 2*pi/length(uu)
end

#= Recursive computation of H3 in spectral space. =#
function ham3direct(uhat::Vector{Complex{Float64}})
	if length(uhat) == 2
		return real( conj(uhat[end]) * uhat[1]^2 )
	else
		H3 = ham3direct(uhat[1:end-1])
		usum = sum( uhat[nn]*uhat[end-nn] for nn=1:lastindex(uhat)-1 )
		H3 += real( conj(uhat[end]) *  usum)
		return H3
	end
end

#---------- High-level Hamiltonian Routines ----------#
# Sampled state
struct RandList
	rvar::Array; H2::Vector; H3::Vector; 
end
# List of constants
struct ConstantList
	C2::Float64; C3::Float64; D0::Float64
end
#= Compute the constanst C2 and C3. =#
function C2C3(eps0::Float64, del0::Float64, lamfac::Int)
	C2 = (2/3) * pi^2 * del0/lamfac^2
	C3 = 1.5 * sqrt(pi) * eps0/del0
	return C2,C3
end
#= Compute the Hamiltonian. H3 and H2 are real and can be numbers or vectors. =#
function hamiltonian(H2, H3, C2::Float64, C3::Float64, D0::Float64)
	return C2*D0^(1/2)*H2 - C3*D0^(-3/2)*H3
end
#= Compute the upstream Hamiltonian.=#
function hamup(rr::RandList, cc::ConstantList)
	return hamiltonian(rr.H2, rr.H3, cc.C2, cc.C3, 1.)
end
#= Compute the downstream Hamiltonian.=#
function hamdn(rr::RandList, cc::ConstantList)
	return hamiltonian(rr.H2, rr.H3, cc.C2, cc.C3, cc.D0)
end
#= Mean of ham under Gibbs measure hamgibbs. =#
function gibbs_mean(quantity, hamgibbs, theta::Float64)
	return dot(exp.(-theta*hamgibbs), quantity) / sum(exp.(-theta*hamgibbs))
end
#= Compute the skewness of the displacement, u or eta. =#
function skewu(H3, hamgibbs, theta::Float64)
	meanH3 =  gibbs_mean(H3, hamgibbs, theta)
	skewu = 3*pi^(1/2)*meanH3
	return skewu
end
#---------------------------------------#

#---------- Sampling Routines ----------#
#= Get uhat from the array of random values. =#
function getuhat(rvar::Array{Float64}, nn::Int)
	uhat = rvar[:,1,nn] + im*rvar[:,2,nn]
	return uhat/sqrt(energy(uhat))
end
#= Sample a single sweep from a uniform distribution on the hypershpere E=1. =#
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
#= Sample several sweeps for better memory usage. =#
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

#---------- Highest-level Routines ----------#
#= Calculate the transfer function:
Given an upstream theta, enforce the matching condition to compute the downstream theta.=#
function transfun(randfile::AbstractString, lamfac::Int)
	# Set Parameters.
	thup = 0.: 0.5 : 30.
	eps0 = 0.017
	del0 = 0.23
	D0 = 0.24
	# Load data.
	randfile = string(data_folder(), randfile)
	rr, nmodes, nsweeps = load(randfile, "rr", "nmodes", "nsweeps")
	# Set Constants.
	C2,C3 = C2C3(eps0,del0,lamfac)
	cc = ConstantList(C2,C3,D0)
	# Compute the list of Hamiltonians using precomputed H2 and H3.
	hdn = Float32.(hamdn(rr,cc))
	hup = Float32.(hamup(rr,cc))
	# Initialize variables.
	nth = length(thup)
	thdn, skup, skdn = [zeros(nth) for idx=1:3]
	guess = thup[1]
	# Loop over the theta values to compute the transfer function.
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
	println("\nCompleted sampling with nmodes = ", nmodes, ", nsweeps = ", nsweeps, " lamfac = ", lamfac)
	println("CPU time = ", sig(cputime/60,3), " minutes.\n\n")	
	savefile = string(data_folder(), "thvars-", string(nmodes), "z-", 
		string(lamfac), "-", string(nsweeps), ".jld")
	save(savefile, "thup", thup, "thdn", thdn, "skup", skup, "skdn", skdn, 
		"nmodes", nmodes, "lamfac", lamfac, "nsweeps", nsweeps, "cputime", cputime)
end

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
