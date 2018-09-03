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






#= Compute the expected value of the downstream Hamiltonian under
either the upstream or downstream Gibbs measure with given theta. =#
function meanham(H3vec::Vector{Float64}, H2vec::Vector{Float64}, 
		theta::Float64, E0::Float64, D0::Float64, gibbsup::Bool)
	ham_dn_mean, norm_const = 0.,0.
	for nn=1:endof(H3vec)
		# Compute the upstream and downstream Hamiltonians.
		hamup = sqrt(E0)*H3vec[nn] - H2vec[nn]
		hamdn = D0^(-13/4)*sqrt(E0)*H3vec[nn] - D0^(3/2)*H2vec[nn]
		# Decide whether to use Gup or Gdn
		gibbsup? ham = hamup : ham = hamdn
		# Compute the mean of Hdn under Gup.
		ham_dn_mean += exp(-theta*ham) * hamdn
		norm_const += exp(-theta*ham)
	end
	return ham_dn_mean/norm_const
end

#= Enforce the statistical matching condition. =#
function matchmean(nmodes::Int, nsamptot::Int, E0::Float64, D0::Float64)
	# Set the upstream inverse temperatures to use.
	thup_vec = -.4:0.05:0.1
	thdn_vec = zeros(Float64,0)
	# Sample H3 and H2 from a microcanonical distribution.
	H3vec, H2vec, rvar = microcan(nmodes,nsamptot)
	# Define a function for the downstream mean of the downstream Hamiltonian.
	meanham_dn(theta_dn::Float64) = meanham(H3vec,H2vec,theta_dn,E0,D0,false)
	# For each theta_up, find the corresponding theta_dn by matching the mean.
	for theta_up in thup_vec
		# Compute the upstream mean of the downstream Hamiltonian.
		mean_up = meanham(H3vec,H2vec,theta_up,E0,D0,true)
		# Find the theta_dn via a root find.
		meandiff(theta_dn::Float64) = meanham_dn(theta_dn) - mean_up
		theta_dn = find_zero(meandiff, theta_up, Order1())
		push!(thdn_vec,theta_dn)
	end
	plot(theta_up,theta_dn,"-.",xlabel="theta_up",ylabel="theta_dn")
	return 
end

matchmean(10, 1*10^5, 4., 0.5)



