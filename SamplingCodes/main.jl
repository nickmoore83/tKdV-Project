include("routines.jl")

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

function HdnGup(H3vec::Vector{Float64}, H2vec::Vector{Float64}, 
		thup::Float64, E0::Float64, D0::Float64)
	Hdn_mean, normconst = 0.,0.
	for nn=1:endof(H3vec)
		# Compute the upstream and downstream Hamiltonians.
		Hup = sqrt(E0)*H3vec[nn] - H2vec[nn]
		Hdn = D0^(-13/4)*sqrt(E0)*H3vec[nn] - D0^(3/2)*H2vec[nn]
		# Compute the mean of Hdn under Gup.
		Hdn_mean += exp(-thup*Hup) * Hdn
		normconst += exp(-thup*Hup)
	end
	return Hdn_mean/normconst
end

function main()
	H3vec, H2vec, rvar = microcan(nmodes,nsamples)
	HdnGup()





