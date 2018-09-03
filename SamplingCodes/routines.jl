using Distributions
using Roots
using Plots

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
	# Compute H3 and H2 for each to estimate the acceptance rate.
	# TO DO: parallelize this for loop!!! It is the bottleneck!
	for nn=1:nsamples
		uhat = getuhat(rvar,nn)
		H3vec[nn] = ham3(uhat)
		H2vec[nn] = ham2(uhat)
		if mod(nn, 10^4) == 0
			println("Initial sampling is ", signif(100*nn/nsamples,3), "% completed.")
		end
	end
	return H3vec, H2vec, rvar
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
#= =#
function gibbs_sample(rvar::Array{Float64}, H3all::Vector{Float64}, H2all::Vector{Float64},
		E0::Float64, D0::Float64, theta::Float64, savemicro::Bool)
	# Set some parameters.
	maxusamples = 2*10^5
	nmodes = size(rvar)[1]
	nsamptot = size(rvar)[3]
	println("\nSampling with D0 = ", signif(D0,2), " and theta = ", signif(theta,2))
	# Determine the acceptance rate based on the Hamiltonian.
	hamvec = D0^(-13/4)*sqrt(E0)*H3all - D0^(3/2)*H2all
	accept_vec = exp.(-theta * hamvec)
	accept_vec *= 1/(maximum(accept_vec))
	# Given the normalized acceptance rate, accept or reject each.
	univar = rand(nsamptot)
	H3acc, H2acc = [zeros(Float64,0) for nn=1:2]
	savemicro? uhacc = zeros(Complex128,nmodes,maxusamples) : 0
	count = 0
	for nn=1:nsamptot
		# Decide to accept or not based on univar and accept_vec
		if univar[nn] <= accept_vec[nn]
			count += 1
			push!(H3acc, H3all[nn]); push!(H2acc, H2all[nn])
			if savemicro && count <= maxusamples
				uhacc[:,count] = getuhat(rvar,nn)
			end
		end
		# Print progress.
		if mod(nn, 10^4) == 0
			println("Acceptance/rejection loop is ", signif(100*nn/nsamptot,3), "% completed.")
		end
	end
	# Remove unused entries of uhacc.
	savemicro? uhacc = uhacc[:,1:min(count,maxusamples)] : uhacc = []

	# Write to data files
end


