

#= Might still finish.
Find indices corresponding to the largest values of H3 
and smallest values of H2.=#
function findind(H3::Vector{Float64},H2::Vector{Float64})
	# FILL IN
	return H3ind, H2ind
end



#----- GOOD ROUTINES THAT WORKED BUT THAT HAVE BEEN SUPERSEDED ----#
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

