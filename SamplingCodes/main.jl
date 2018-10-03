include("routines.jl")

#---------- Reading and writing routines ----------#
# Extract the parameters
function extractparams(params::Vector)
	nmodes, nsamp, nsweeps = Int(params[1]), Int(params[2]), Int(params[3])
	amp, D0 = params[4:5]
	lamfac = Int(params[6]) 
	thmin, thmax, dth = params[7:9]
	thup_vec = collect(thmin:dth:thmax)
	nthetas = endof(thup_vec)
	return nmodes, nsamp, nsweeps, amp, D0, lamfac, thup_vec, nthetas
end
# Write the macrostate data.
function write_mac_data(accstate::AcceptedState, theta::Float64, 
		suffix::AbstractString, run_number::Int)
	# Preliminaries.
	foldername = datafolder(run_number)
	macfile = string(foldername,"mac",suffix,".txt")
	nacc = accstate.naccepted
	# Get coarse information from microstates.
	micmax, macmax = maxparams()
	uhat = accstate.uhat[:,1:min(nacc,micmax)]
	uhavg = getuhavg(uhat)
	# Save the macrostates.
	println("Saving the macrostates ", suffix)
	label1 = "# Macrostate data"
	label2 = "# Basic information: theta, number of accepted samples"
	label3 = "# Computed data: vectors of accepted H3 and H2, vector of mean uhat per mode"
	macdata = [label1; label2; theta; nacc;
				label3; accstate.H3[1:nacc]; accstate.H2[1:nacc]; uhavg]
	writedata(macdata, macfile)
	# Save the microstates if desired.
	# Note: this step is expensive because it requires transforming to physical space.
	if savemicro()
		println("Saving the microstates ", suffix)
		uacc = getuacc(uhat)
		micfile = string(foldername,"mic",suffix,".txt")
		writedata(uacc, micfile)
	end
end
# Write all the data, including basic, microstate, and mincrostate.
function write_all_data(params::Vector, thdn_vec::Vector{Float64},
		accstate::Array{AcceptedState}, cputimes::Vector{Float64}, run_number::Int)
	foldername = datafolder(run_number)
	# Write the basic data.
	basicfile = string(foldername,"basic.txt")
	nmodes, nsamp, nsweeps, amp, D0, lamfac, thup_vec, nthetas = extractparams(params)
	label1 = "# Input parameters: "
	label2 = "# CPU times for matching mean and for sampling (mins)"
	label3 = "# Inverse temperature data: number of thetas, theta_ups and theta_dns"
	basicdata = [label1; params; label2; cputimes;
		label3; nthetas; thup_vec; thdn_vec]
	writedata(basicdata,basicfile)
	# Write the macro data.
	for nn=1:nthetas
		write_mac_data(accstate[nn,1], thup_vec[nn], string("up",nn), run_number)
		write_mac_data(accstate[nn,2], thdn_vec[nn], string("dn",nn), run_number)
	end
end
#---------------------------------------#

#---------- Routines to enforce matching condition ----------#
#= Compute the expected value of the downstream Hamiltonian under
either the upstream or downstream Gibbs measure with given theta. =#
function meanham(H3vec::Vector{Float64}, H2vec::Vector{Float64}, 
		amp::Float64, D0::Float64, lamfac::Int, theta::Float64, gibbsup::Bool)
	hamdn = getham(H3vec,H2vec,amp,D0,lamfac)
	gibbsup? ham=getham(H3vec,H2vec,amp,1.,lamfac) : ham = hamdn
	return dot(exp.(-theta*ham), hamdn) / sum(exp.(-theta*ham))
end
#= Determine the downstream thetas that satisfy the statistical matching condition. =#
function matchmean(nmodes::Int, nsamp::Int, 
		amp::Float64, D0::Float64, lamfac::Int, thup_vec::Vector{Float64})
	# Preliminaries
	nthetas = endof(thup_vec)
	thdn_vec = zeros(Float64,nthetas)
	# Sample H3 and H2 from a microcanonical distribution.
	rset = microcan(nmodes,nsamp)
	# For each thup, find the corresponding thdn by matching the mean.
	meanham_dn(theta_dn::Float64) = meanham(rset.H3,rset.H2,amp,D0,lamfac,theta_dn,false)
	for nn = 1:nthetas
		# Determine thdn by finding a root of the difference of means.
		mean_up = meanham(rset.H3,rset.H2,amp,D0,lamfac,thup_vec[nn],true)
		meandiff(theta_dn::Float64) = meanham_dn(theta_dn) - mean_up
		thdn_vec[nn] = find_zero(meandiff, thup_vec[nn], Order1())
	end
	return thdn_vec, rset
end
#---------------------------------------#

#= Main routine to enforce the statistical matching condition. =#
function main(run_number::Int=0)
	# Preliminaries.
	newfolder(datafolder(run_number))
	paramsfile = string("params",run_number,".txt")
	params = readvec(paramsfile)
	nmodes, nsamp, nsweeps, amp, D0, lamfac, thup_vec, nthetas = extractparams(params)
	# Determine thdn to match the means.
	println("Enforcing the statistical matching condition.")
	cput_match = @elapsed (thdn_vec, rset) = matchmean(nmodes,nsamp,amp,D0,lamfac,thup_vec)
	cput_match = signif(cput_match/60,2)
	println("The CPU time for enforcing matching condition is ", cput_match, " minutes.")
	#plt = plot(thup_vec,thdn_vec, xlabel="theta_up",ylabel="theta_dn"); display(plt)	
	# Initialize the accepted set of states.
	accstate = Array{AcceptedState}(nthetas,2)
	for nn=1:nthetas
		accstate[nn,1] = new_acc_state(nmodes)
		accstate[nn,2] = new_acc_state(nmodes)
	end
	#= Define a function to sample from the Gibbs distributions
	for all values of upstream and downstream values of theta. =#
	function gibbs_sample_updn(rset::RandSet, accstate::Array{AcceptedState})
		for nn = 1:nthetas
			gibbs_sample!(rset, accstate[nn,1], amp,1.,lamfac,thup_vec[nn])
			gibbs_sample!(rset, accstate[nn,2], amp,D0,lamfac,thdn_vec[nn])
		end
	end
	tm0 = time()
	# Sample using rvar, H3, and H2 from the matchmean computation.
		println("\n\n\nBeginning sweep 1 of sampling.")
	gibbs_sample_updn(rset,accstate)
	# Take several additional sweeps sampling from the Gibbs distributions.
	for sweep = 2:nsweeps
		println("\n\n\nBeginning sweep ", sweep, " of sampling.")
		rset = microcan(nmodes,nsamp)
		gibbs_sample_updn(rset,accstate)
	end
	cput_sample = signif((time()-tm0)/60, 2)
	println("\n\nCompleted the Gibbs sampling phase.")
	println("The CPU time for sampling is ", cput_sample, " minutes.")
	cputimes = [cput_match, cput_sample]
	write_all_data(params,thdn_vec,accstate,cputimes,run_number)
end

