
#----- GOOD ROUTINES THAT WORKED BUT THAT HAVE BEEN SUPERSEDED ----#

#= Compute the two Hamiltonians (upstream and downstream)
as sampled from the Gibbs0 distribution. =# 
function sample_gibbs(nmodes::Int, nsamples::Int, invtemp::Float64, 
		E0::Float64, D0::Float64, suffix::AbstractString)
	savemicro = true
	println("\n\n")
	# Sample from a Gibbs distribution.
	cputime = @elapsed (H3vec,H2vec,uhacc) = gibbs(nmodes,nsamples,invtemp,E0,D0,savemicro)
	cputime = signif(cputime/60,2)	
	naccept = endof(H3vec)
	acceptrate = signif(100*naccept/nsamples,2)
	# Convert each accepted uhat to physical space for analysis.
	uacc = getuacc(uhacc)
	uhavg = getuhavg(uhacc)
	# Write Hamiltonian data to a file.
	println("\nWriting data to output file.")
	hamfile = string("../SamplingData/ham",suffix,".txt")
	label1 = "# Input parameters: nmodes, nsamples, invtemp, E0, D0"
	label2 = "# Calculated parameters: number of accepted samples, acceptance rate (%), CPU time (minutes)"
	label3 = "# Computed data: H3, H2, uhavg"	
	hamdata = [label1; nmodes; nsamples; invtemp; E0; D0;
				label2; naccept; acceptrate; cputime; label3; H3vec; H2vec; uhavg]	
	writedata(hamdata, hamfile)
	# Write uacc data to a file.
	ufile = string("../SamplingData/u",suffix,".txt")
	writedata(uacc, ufile)
	# Print the final CPU time.
	println("The CPU time is ", signif(cputime,2), " minutes.")
	return
end

# Testing case:
#sample_gibbs(10, 1*10^5, -0.05, 4., 1.0, "up")
# Three main runs:
#sample_gibbs(10, 1*10^6, -0.0, 4., 1.0, "up")
#sample_gibbs(10, 2*10^7, -0.3, 4., 1.0, "up")
#sample_gibbs(10, 2*10^7, -0.3, 4., 0.5, "dn")

#= Sample from a Gibbs distribution with arbitrary inverse temperature. =#
function gibbs(nmodes::Int, nsamples::Int, invtemp::Float64=0., 
		E0::Float64=1., D0::Float64=1., savemicro::Bool=false)
	# Sample H3 and H2 from the microcanonical distribution.
	H3all,H2all,rvar = microcan(nmodes,nsamples)
	# Determine the accpetance rate based on the Hamiltonian.
	aratevec = zeros(Float64,nsamples)
	maxaccept = 0.
	for nn=1:nsamples
		# Compute the Hamiltonian (for upstream D0=1)
		ham = D0^(-13/4)*sqrt(E0)*H3all[nn] - D0^(3/2)*H2all[nn]
		# Compute the unnormalized acceptance rate and find the maximum.
		aratevec[nn] = exp(-invtemp * ham)
		maxaccept = max(maxaccept, aratevec[nn])
	end
	# Now that the normalization is known, accept or reject each.
	H3acc,H2acc = [zeros(Float64,0) for nn=1:2]
	savemicro? uhacc = zeros(Complex128,nmodes,nsamples):0
	count = 0
	for nn=1:nsamples
		# Accept or reject based on a uniform random variable.
		acceptrate = aratevec[nn]/maxaccept
		univar = rand()
		if univar <= acceptrate
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
		# Print progress.
		if mod(nn, 10^4) == 0
			println("Acceptance/rejection loop is ", signif(100*nn/nsamples,3), "% completed.")
		end
	end
	# Remove unused entries of uhacc.
	savemicro? uhacc = uhacc[:,1:count] : uhacc = []
	# Print the overall acceptance rate.
	println("The overall acceptance rate was ", signif(100*count/nsamples,2), "%.")
	return H3acc, H2acc, uhacc
end

# OBSELETE since I wrote gibbs for arbitrary invers temperature.
#= Sample from a microcanonical distribution with fixed energy E=1.0,
i.e. a Gibbs distribution with zero inverse temperature. =#
function gibbs0(nmodes::Int, nsamples::Int)
	# Get all the random samples at once.
	rvar = randn(nmodes,2,nsamples)
	H3 = zeros(Float64,nsamples)
	H2 = zeros(Float64,nsamples)
	# Compute H3 and H2 for each microstate. This is the bottleneck.
	# TO DO: parallelize this for loop!!!
	for nn=1:nsamples
		uhat = getuhat(rvar,nn)
		H3[nn] = ham3(uhat)
		H2[nn] = ham2(uhat)
		if mod(nn, 10^5) == 0
			println("Computation is ", signif(100*nn/nsamples,3), "% completed.")
		end
	end
	return H3, H2
end
# OBSELETE
#= Compute H3 using direction summation in spectral space. =#
function ham3sum(uhat::Vector{Complex128})
	nmodes = endof(uhat)
	H3 = 0.
	# Compute H3 directly with a double sum over k1 and k2.
	for k1 = 1:nmodes
		for k2 = 1:nmodes
			#= Case: k1 and k2 have same sign. 
			Consider k1,k2>0, then add CC via 2*Real =#
			k3 = k1 + k2
			if k3 <= nmodes
				H3 += 2*real(uhat[k1]*uhat[k2]*conj(uhat[k3]))
			end
			#= Case: k1 and k2 have opposite signs.
			Consider k1>0 & k2<0, then add CC via 2*Real =#
			k3 = k2-k1
			if (abs(k3) <= nmodes)
				# Case abs(k2)>abs(k1)
				if k3 > 0
					H3 += 2*real(uhat[k1]*conj(uhat[k2])*uhat[k3])
				# Case abs(k1)>abs(k2)
				elseif k3 < 0
					H3 += 2*real(uhat[k1]*conj(uhat[k2])*conj(uhat[-k3]))
				end
			end
		end
	end
	# Normalize and return as output.
	return H3/6
end
# OBSELETE
#= Sample from a microcanonical distribution with fixed energy E=1.0,
i.e. a Gibbs distribution with zero inverse temperature. =#
function gibbs0OLD(nmodes::Int)
	rvar = randn(nmodes,2)
	uhat = rvar[:,1] + im*rvar[:,2]
	return uhat/sqrt(energy(uhat))
end
#------------------------------------------------------------------#

#----- IDEAS THAT TURNED OUT TO BE BAD ----#
# INITIAL IDEA FOR GIBBS MEASURE; PROBABLY WILL NOT USE.
#= Sample from a Gibbs distribution with arbitrary inverse temperature. 
This algorithm was for taking a single sample. =#
function gibbs(nmodes::Int, invtemp::Float64=0.)
	accepted = false
	while(accepted = false)
		# Sample uhat from Gibbs0 and compute H3 and H2.
		uhat = gibbs0(nmodes)
		H3 = ham3(uhat)
		H2 = ham2(uhat)
		# Set the coefficients and compute the Hamiltonian
		#A3 = sqrt(E0)*lambda^(-3/2)
		#A2 = lambda^(-3)
		A3 = 1.
		A2 = 1.
		ham = A3*H3 - A2*H2
		# Compute the accptance probability using the inverse temperature.
		acceptprob = exp(-invtemp * ham)
		# NEED TO NORMALIZE....
		uvar = rand()
		if uvar <= acceptprob
			accepted = true
		end
	end
end
# PROBABLY WILL NOT WORK OUT
#= This was a good idea, but I think it will not work after all because 
I would have to figure out how the distributions of H3 and H2 depend on
the number of modes. Dang! That is the only downside. =#
function isexceptional(H3::Float64,H2::Float64,nmodes::Int)
	H3verybig = 2.
	H3big = 1.
	H2verysmall = 1.
	H2small = 2.
	# If either H3 or H2 is very extreme, return true.
	if ((H3 > H3verybig) | (H2 < H2verysmall))
		return true
	# If both H3 and H2 are pretty extreme, return true.
	elseif (H3 > H3big & H2 < H2small)
		return true
	else
		return false
	end
end
#------------------------------------------------------------------#

#----- MISC ----#
		# Print progress.
		if mod(nn, 10^4) == 0
			println("Acceptance/rejection loop is ", signif(100*nn/nsamptot,3), "% completed.")
		end

function matchmean(nmodes::Int, nsamptot::Int, E0::Float64, D0::Float64)

