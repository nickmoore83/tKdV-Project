using Distributions
using Roots
using Plots

#---------- Data Types ----------#
# Fix some parameters.
function maxparams()
	micmax = 2*10^5
	macmax = 5*10^6
	return micmax, macmax
end
# Sampled state
type RandSet
	H3::Vector{Float64}; H2::Vector{Float64}; rvar::Array{Float64}
end
# Accepted state
type AcceptedState
	H3::Vector{Float64}; H2::Vector{Float64}; uhat::Array{Complex128}; naccepted::Int
end
# Initiate a new StateAccepted
function new_acc_state(nmodes::Int)
	micmax, macmax = maxparams()
	return AcceptedState(zeros(Float64,macmax), zeros(Float64,macmax), 
		zeros(Complex128,nmodes,macmax), 0)
end
#---------------------------------#

#---------- IO Routines ----------#
function readvec(file::AbstractString)
	iostream = open(file, "r")
	vec1 = readdlm(iostream)[:,1]
	close(iostream)
	return vec1
end
function writedata(data::Array, filename::AbstractString)
	iostream = open(filename, "w")
	writedlm(iostream, data)
	close(iostream)
end
function datafolder()
	return "../SamplingData/run/"
end
function newfolder(foldername::AbstractString)
	isdir(foldername)? rm(foldername; recursive=true) : 0
	mkdir(foldername)
end
#---------------------------------#

#---------- Real FFT Routines ----------#
#= Compute the FFT between physical and spectral space assuming
1) the signal uu is real and 2) uu has zero mean (i.e. momentum M=0)
Note: uu has length npoints = 2*nmodes. =#
#= Basic realfft to go from uu to uhat. =#
function realfft(uu::Vector{Float64})
	uhat = rfft(uu)/endof(uu)
	assert(abs(uhat[1])/max.(abs,uhat) < 1e-6) 
	return uhat[2:end]
end
#= Basic irealfft to go from uhat to uu (no upsampling). =#
function irealfft(uhat::Vector{Complex128})
	nmodes = endof(uhat)
	uu = irfft([0; uhat], npoints(nmodes))
	return uu*endof(uu)
end
#= Upsampled version of irealfft. =#
function ifftup(uhat::Vector{Complex128})
	uhat = [uhat; zeros(eltype(uhat), endof(uhat))]
	return irealfft(uhat)
end
#= Valid choices for npts are 2*nmodes or 2*nmodes+1, 
but I need a function to keep my usage consistent.
Note: for H3test npts = 2*modes is slightly more accurate. =#
function npoints(nmodes::Int)
	return 2*nmodes
end
#---------------------------------------#

#---------- Hamiltonian Routines ----------#
#= In all routines, the physical signal is assumed to be real 
with zero momentum, M = int u dx = 0. Therefore, the Fourier transform 
only requires the modes k=1,...,nmodes. The negative modes are given 
by the complex conjugates (CC), and the k=0 mode vanishes due to M=0. 
I used the shorthand convention int = (1/2*pi) int_0^{2*pi}. 
and nmodes is the same as Lambda. =#

#= Compute the energy, E = 1/2 int u^2 dx. =#
function energy(uhat::Vector{Complex128})
	return norm(uhat)^2
end
# Compute H2 = 1/2 int u_x^2 dx
function ham2(uhat::Vector{Complex128})
	kvec = 1:endof(uhat)
	uxhat = im*kvec.*uhat
	return energy(uxhat)
end
#= Compute H3 using FFT to physical space, H3 = 1/6 int u^3 dx. =#
function ham3(uhat::Vector{Complex128})
	uu = ifftup(uhat)
	return sum(uu.^3) / (6*endof(uu))
end
#---------------------------------------#

#---------- Short Sampling Routines ----------#
#= Get uhat from the array of random values. =#
function getuhat(rvar::Array{Float64}, nn::Int)
	uhat = rvar[:,1,nn]+im*rvar[:,2,nn]
	return uhat/sqrt(energy(uhat))
end
#= Sample from a microcanonical distribution (Gibbs with zero inverse temperature). =#
function microcan(nmodes::Int, nsamples::Int)
	# Get all the random samples at once.
	rvar = randn(nmodes,2,nsamples)
	# Allocate space.
	H3vec, H2vec = [zeros(Float64,nsamples) for nn=1:2]
	println("\nSampling from microcanonical distribution.")
	# Compute H3 and H2 for each to estimate the acceptance rate.
	# TO DO: parallelize this for loop!!! It is the bottleneck!
	for nn=1:nsamples
		uhat = getuhat(rvar,nn)
		H3vec[nn] = ham3(uhat)
		H2vec[nn] = ham2(uhat)
		if mod(nn, 10^4) == 0
			println("Microcanonical sampling is ", signif(100*nn/nsamples,3), "% completed.")
		end
	end
	return RandSet(H3vec,H2vec,rvar)
end
#= Convert each accepted uhat to physical space for analysis. =#
function getuacc(uhacc::Array{Complex128})
	nmodes = size(uhacc)[1]
	nsamples = size(uhacc)[2]
	npts = npoints(nmodes)
	uacc = zeros(Float64,npts,nsamples)
	for nn=1:nsamples
		uacc[:,nn] = irealfft(uhacc[:,nn])
	end
	return uacc
end
#= Compute the ensemble average of abs(uhat) for each mode =#
function getuhavg(uhacc::Array{Complex128})
	nmodes = size(uhacc)[1]
	uhavg = zeros(Float64,nmodes)
	for kk=1:nmodes
		uhavg[kk] = mean(abs.(uhacc[kk,:]))
	end
	return uhavg
end
#---------------------------------------#

#---------- Main Sampling Routines ----------#
#= Sample from a Gibbs distribution with non-zero theta, 
given that H3 and H2 have already been sampled from microcanonical distribution. =#
function gibbs_sample!(rset::RandSet, accstate::AcceptedState, 
		E0::Float64, D0::Float64, theta::Float64)
	# If macmax is already exceeded, then return immediately.	
	micmax, macmax = maxparams()
	if accstate.naccepted >= macmax; return; end
	# Preliminaries	
	H3all = rset.H3[:]
	H2all = rset.H2[:]
	rvar = rset.rvar[:,:,:]
	nmodes = size(rvar)[1]
	nsamp = size(rvar)[3]
	# Determine the acceptance rate based on the Hamiltonian.
	println("\nSampling from Gibbs distribution with D0 = ", 
		signif(D0,2), " and theta = ", signif(theta,2))
	hamvec = D0^(-13/4)*sqrt(E0)*H3all - D0^(3/2)*H2all
	accept_vec = exp.(-theta * hamvec)
	accept_vec *= 1/(maximum(accept_vec))
	# With the normalized acceptance rate, decide to accept/reject each.
	univar = rand(nsamp)
	for nn=1:nsamp
		accstate.naccepted >= macmax? break : 0
		if univar[nn] <= accept_vec[nn]
			accstate.naccepted += 1
			mm = accstate.naccepted
			accstate.H3[mm] = H3all[nn]
			accstate.H2[mm] = H2all[nn]
			if mm <= micmax
				accstate.uhat[:,mm] = getuhat(rvar,nn)
			end
		end
	end
	println("Completed acceptance/rejection phase.")
	return
end
