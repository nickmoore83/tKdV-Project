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

#---------- Real FFT Routines ----------#
#= Compute the FFT between physical and spectral space assuming
1) the signal uu is real and 2) uu has zero mean (i.e. momentum M=0)
Note: uu has length npoints = 2*nmodes or 2*nmodes + 1. =#
#= Valid choices for npoints are 2*nmodes or 2*nmodes+1, 
but I need a function to keep my usage consistent.
Note: for H3test npts = 2*modes is slightly more accurate. =#
function npoints(nmodes::Int)
	return 2*nmodes
end
#= Basic realfft to go from uu to uhat. =#
function realfft(uu::Vector{Float64})
	uhat = rfft(uu)/length(uu)
	@assert(abs(uhat[1])/maximum(abs,uhat) < 1e-6) 
	return uhat[2:end]
end
#= Basic irealfft to go from uhat to uu (no upsampling). =#
function irealfft(uhat::Vector{Complex{Float64}})
	nmodes = length(uhat)
	uu = irfft([0; uhat], npoints(nmodes))
	return uu*length(uu)
end
#= Upsampled version of irealfft. =#
function ifftup(uhat::Vector{Complex{Float64}})
	uhat = [uhat; zeros(eltype(uhat), length(uhat))]
	return irealfft(uhat)
end
#---------------------------------------#

#---------- Hamiltonian Routines ----------#
#= In all routines, the physical signal is assumed to be real 
with zero momentum, M = int u dx = 0. Therefore, the Fourier transform 
only requires the modes k=1,...,nmodes. The negative modes are given 
by the complex conjugates (CC), and the k=0 mode vanishes due to M=0. =#
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
function hamiltonian(H2, H3, C2, C3, D0)
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
function meanham(ham, hamgibbs, theta)
	return dot(exp.(-theta*hamgibbs), ham) / sum(exp.(-theta*hamgibbs))
end
#= Compute the skewness of the displacement, u or eta. =#
function skewu(H3, hamgibbs, theta)
	meanH3 =  dot(exp.(-theta*hamgibbs), H3) / sum(exp.(-theta*hamgibbs))
	skewu = 3*pi^(1/2)*meanH3
	return skewu
end
#---------------------------------------#

#---------- Sampling Routines ----------#
#= Get uhat from the array of random values. =#
function getuhat(rvar::Array{Float64}, nn::Int)
	uhat = rvar[:,1,nn]+im*rvar[:,2,nn]
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
		H2vec[nn] = ham2(uhat); H3vec[nn] = ham3(uhat)
	end
	rlist = RandList(Float32.(rvar), H2vec, H3vec)
	return rlist
end
#= Sample several sweeps from a uniform distribution on the hypershpere E=1.
Uses several sweeps manages the memory better. 
Note: zerolast zeroes the last mode as is done in the Matlab DNS. 
Note: savemicro not yet addressed at all. =#
function uniform_sample(nmodes::Int, nsweeps::Int, 
		zerolast::Bool=true, savemicro::Bool=false)
	samp_per = 10^5
	totsamp = nsweeps*samp_per
	H2vec, H3vec = [zeros(Float32,totsamp) for nn=1:2]
	for nn = 1:nsweeps
		rlist = sample_one_sweep(nmodes, samp_per, zerolast)
		indrange = (nn-1)*samp_per+1 : nn*samp_per
		H2vec[indrange] = Float32.(rlist.H2)
		H3vec[indrange] = Float32.(rlist.H3)
		println("Sweep ",nn,", ",sig(100*nn/nsweeps,3),"% completed.")
	end
	zerolast ? zstr = "z" : zstr = ""
	savefile = string("rand-",string(nmodes),zstr,"-",string(nsweeps),".jld")
	save(savefile, "rr", RandList([],H2vec,H3vec),
		"nmodes", nmodes, "nsweeps", nsweeps, "totsamp", totsamp)

	# OLD RELEVANT CODE SNIPPET FOR SAVEMICRO:
	#savemicro ? rsave = Float32.(rvar) : rsave = []
end
#---------------------------------------#





function writedata(data::Array, filename::AbstractString)
	iostream = open(filename, "w")
	writedlm(iostream, data)
	close(iostream)
end

function readvec(file::AbstractString)
	iostream = open(file, "r")
	vec1 = readdlm(iostream, comments=true, comment_char='#')[:,1]
	close(iostream)
	return vec1
end

# Significant digits
function sig(var, sigdig)
	return round(var,sigdigits=sigdig)
end