include("routines.jl")

#= Compute the two Hamiltonians (upstream and downstream)
as sampled from the Gibbs0 distribution. =# 
function main(nmodes::Int, nsamples::Int, invtemp::Float64, E0::Float64, D0::Float64)
	savemicro = true
	println("\n\n")
	# Sample from a Gibbs distribution.
	cputime = @elapsed( (H3up,H2up,uhup,H3dn,H2dn,uhdn) = 
				gibbs(nmodes,nsamples,invtemp,E0,D0,savemicro) )
	cputime = signif(cputime/60,2)

	# Function to save the data.
	function save_ham_data(H3::Vector{Float64}, H2::Vector{Float64}, 
			uhat::Vector{Float64}, suffix::AbstractString)
		# Calculate the acceptance rate.
		naccept = endof(H3)
		acceptrate = signif(100*naccept/nsamples,2)
		# Convert each accepted uhat to physical space for analysis.
		uu = getuacc(uhat)
		uhavg = getuhavg(uhat)
		# Write Hamiltonian data to a file.
		println("\nWriting data to output file.")
		label0 = string("# ", suffix, "stream data")
		label1 = "# Input parameters: nmodes, nsamples, invtemp, E0, D0"
		label2 = "# Calculated parameters: number of accepted samples, acceptance rate (%), CPU time (minutes)"
		label3 = "# Computed data: H3, H2, uhavg"
		hamdata = [label0; label1; nmodes; nsamples; invtemp; E0; D0;
				label2; naccept; acceptrate; cputime;
				label3; H3vec; H2vec; uhavg]
		hamfile = string("../SamplingData/hamdata",suffix,".txt")
		writedata(hamdata, hamfile)
		# Write uacc data to a file.
		#udata = [uacc [real(uhacc); imag(uhacc)] ]
		ufile = string("../SamplingData/udata",suffix,".txt")
		writedata(uu, ufile)
	end
	# Save the data.
	save_ham_data(H3up,H2up,uhup,"up")
	save_ham_data(H3dn,H2dn,uhdn,"dn")
	# Print the final CPU time.
	println("The CPU time is ", signif(cputime,2), " minutes.")
	return
end

main(10, 1*10^5, -0.0, 4., 0.5)

#main(10, 1*10^6, -0.0, 4., 0.5)
#main(10, 2*10^7, -0.3, 4., 0.5)
