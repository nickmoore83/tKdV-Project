using Distributions

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
	return
end
#---------------------------------#

#---------- Real FFT Routines ----------#
#= Compute the FFT between physical and spectral space assuming
1) the signal uu is real and 2) uu has zero mean (i.e. momentum M=0)
Note: uu has length npoints = 2*nmodes. =#
#= Basic realfft to go from uu to uhat. =#
function realfft(uu::Vector{Float64})
	uhat = rfft(uu)/endof(uu)
	assert(abs(uhat[1])/maximum(abs,uhat) < 1e-6) 
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

#---------- Sampling Routines ----------#
#= Get uhat from the array of random values. =#
function getuhat(rvar::Array{Float64}, nn::Int)
	uhat = rvar[:,1,nn]+im*rvar[:,2,nn]
	return uhat/sqrt(energy(uhat))
end

#= Sample from a Gibbs distribution with arbitrary inverse temperature. =#
function gibbs(nmodes::Int, nsamples::Int, invtemp::Float64=0., 
		E0::Float64=1., D0::Float64=1., savemicro::Bool=false)
	# Get all the random samples at once.
	rvar = randn(nmodes,2,nsamples)
	# Allocate space.
	H3all, H2all = [zeros(Float64,nsamples) for nn=1:2]
	aupvec, adnvec = [zeros(Float64,nsamples) for nn=1:2]
	maxup, maxdn = 0.,0.
	# Compute H3 and H2 for each to estimate the acceptance rate.
	# TO DO: parallelize this for loop!!! It is the bottleneck!
	for nn=1:nsamples
		uhat = getuhat(rvar,nn)
		H3 = ham3(uhat)
		H2 = ham2(uhat)
		# Compute the Hamiltonian (for upstream D0=1)
		hamup = sqrt(E0)*H3 - H2
		hamdn = D0^(-13/4)*sqrt(E0)*H3 - D0^(3/2)*H2
		# Compute the unnormalized acceptance rate and maximize it.
		accept_up = exp(-invtemp * hamup)
		accept_dn = exp(-invtemp * hamdn)
		maxup = max(maxup, accept_up)
		maxdn = max(maxdn, accept_dn)
		# Save H3 and H2 so as to not recompute.
		H3all[nn] = H3
		H2all[nn] = H2
		aupvec[nn] = accept_up
		adnvec[nn] = accept_dn
		# Print progress.
		if mod(nn, 10^4) == 0
			println("Initial sampling is ", signif(100*nn/nsamples,3), "% completed.")
		end
	end
	# Now that the normalization is known, accept or reject each.
	H3accup, H2accup, H3accdn, H2accdn = [zeros(Float64,0) for nn=1:4]

	savemicro? uhacc = zeros(Complex128,nmodes,nsamples):0 # CHECK

	countup = [0]
	countdn = [0]

	# Little function to accept or reject each sample.
	function accept!(accept_rate::Float64, nn::Int, count::Array{Int}, 
			H3acc::Vector{Float64}, H2acc::Vector{Float64}, uhacc::Vector{Float64})
		univar = rand()
		if univar <= accept_rate
			count += 1
			# Save H3 and H2
			push!(H3acc, H3all[nn])
			push!(H2acc, H2all[nn])
			# Also save the microstate if requested.
			if savemicro
				uhat = getuhat(rvar,nn)
				uhacc[:,count] = uhat[:]
			end
		end	

	for nn=1:nsamples
		# The normalized acceptance rates.
		accept_up = aupvec[nn]/maxup
		accept_dn = adnvec[nn]/maxdn
		# Accept or reject each sample.
		accept!(accept_up, nn, countup, H3accup, H2accup, uhaccup)
		accept!(accept_dn, nn, countdn, H3accdn, H2accdn, uhaccdn)

		# Accept or reject based on a uniform random variable.
		univarup = rand()
		univardn = rand()

		if univar <= accept_up
			countup += 1
			# Save H3 and H2
			push!(H3accup, H3all[nn])
			push!(H2accup, H2all[nn])
			# Also save the microstate if requested.
			if savemicro
				uhat = getuhat(rvar,nn)
				uhacc[:,count] = uhat[:]
			end
		end
		# Print progress.
		if mod(nn, 10^4) == 0
			println("Acceptance/rejection loop is ", signif(100*nn/nsamples,3), "% completed.")
		end
	end
	# Remove unused entries of uhacc.
	savemicro? uhacc = uhacc[:,1:count] : uhacc = []
	# State the overall acceptance rate.
	println("The overall acceptance rate was ", signif(100*count/nsamples,2), "%.")
	return H3acc, H2acc, uhacc
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

