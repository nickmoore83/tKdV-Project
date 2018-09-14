include("routines.jl")

#= Write the basic data to a file. =#
function writebasicdata(params::Vector, nsamptot::Int, 
		cput_match::Float64, cput_sample::Float64, 
		thup_vec::Vector{Float64}, thdn_vec::Vector{Float64})
	foldername = datafolder()
	file = string(foldername,"basic.txt")
	label1 = "# Input parameters: "
	label2 = "# Calculated parameters: nsamptot, CPU time for matching mean and for sampling (mins)"
	label3 = "# Inverse temperature data: number of thetas, theta_ups and theta_dns"
	data = [label1; params; label2; nsamptot; cput_match; cput_sample;
		label3; endof(thup_vec); thup_vec; thdn_vec]
	writedata(data,file)
end

#= Compute the expected value of the downstream Hamiltonian under
either the upstream or downstream Gibbs measure with given theta. =#
function meanham(H3vec::Vector{Float64}, H2vec::Vector{Float64}, 
		E0::Float64, D0::Float64, theta::Float64, gibbsup::Bool)
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

#= Determine the downstream thetas that satisfy the statistical matching condition. =#
function matchmean(nmodes::Int, nsamp::Int, E0::Float64, D0::Float64, thup_vec::Vector{Float64})
	# Preliminaries
	nthetas = endof(thup_vec)
	thdn_vec = zeros(Float64,nthetas)
	# Sample H3 and H2 from a microcanonical distribution.
	H3vec, H2vec, rvar = microcan(nmodes,nsamp)
	# For each thup, find the corresponding thdn by matching the mean.
	meanham_dn(theta_dn::Float64) = meanham(H3vec,H2vec,E0,D0,theta_dn,false)
	for nn = 1:nthetas
		# Determine thdn by finding a root of the difference of means.
		mean_up = meanham(H3vec,H2vec,E0,D0,thup_vec[nn],true)
		meandiff(theta_dn::Float64) = meanham_dn(theta_dn) - mean_up
		thdn_vec[nn] = find_zero(meandiff, thup_vec[nn], Order1())
	end
	return thdn_vec, H3vec, H2vec, rvar
end

#= Main routine to enforce the statistical matching condition. =#
function main(paramsfile::AbstractString="params.txt")
	# Preliminaries.
	savemicro = true
	newfolder(datafolder())
	params = readvec(paramsfile)
	nmodes, nsamp, npasses = Int(params[1]), Int(params[2]), Int(params[3])
	E0, D0, thmin, thmax, dth = params[4:8]
	thup_vec = collect(thmin:dth:thmax)
	# Determine thdn to match the means.
	println("Enforcing the statistical matching condition.")
	cput_match = @elapsed 
		(thdn_vec, H3vec, H2vec, rvar = matchmean(nmodes,nsamp,E0,D0,thup_vec))
	cput_match = signif(cput_match/60,2)
	println("CPU time for enforcing matching condition is ", cputime, " minutes.")
	#plt = plot(thup_vec,thdn_vec, xlabel="theta_up",ylabel="theta_dn"); display(plt)

	#= Define a function to sample from the Gibbs distributions
	for all values of upstream and downstream values of theta. =#
	function gibbs_sample_updn(H3vec::Vector{Float64}, H2vec::Vector{Float64}, rvar::Array{Float64})
		for nn = 1:nthetas
			upsuffix = string("up",nn)
			dnsuffix = string("dn",nn)
			gibbs_sample(H3vec,H2vec,rvar, E0,1.,thup_vec[nn], savemicro,upsuffix)
			gibbs_sample(H3vec,H2vec,rvar, E0,D0,thdn_vec[nn], savemicro,dnsuffix)
		end
	end
	tm0 = time()
	# Sample using rvar, H3, and H2 from the matchmean computation.
	gibbs_sample_updn(H3vec, H2vec, rvar)
	# Take several additional passes sampling from the Gibbs distributions.
	for pass = 2:npasses
		H3vec, H2vec, rvar = microcan(nmodes,nsamp)
		gibbs_sample_updn(H3vec, H2vec, rvar)
	end
	cput_sample = signif((time()-tm0)/60, 2)
	nsamptot = sampper*npasses
	println("The total CPU time is ", cputime, " minutes.")
	writebasicdata(params,nsamptot,cput_match,cput_sample,thup_vec,thdn_vec)
end

