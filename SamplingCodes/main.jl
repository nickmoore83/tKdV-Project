include("routines.jl")

#= Write the basic data to a file. =#
function writebasicdata(params::Vector, nsamptot::Int, cputime::Float64, 
		thup_vec::Vector{Float64}, thdn_vec::Vector{Float64})
	foldername = datafolder()
	file = string(foldername,"basic.txt")
	label1 = "# Input parameters: nmodes, nsampper, npasses, E0, D0, thmin, thmax, dth"
	label2 = "# Calculated parameters: nsamptot, CPU time (mins)"
	label3 = "# Inverse temperature data: number of thetas, theta_ups and theta_dns"
	data = [label1; params; label2; nsamptot; cputime; 
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

#= Main routine to enforce the statistical matching condition. =#
function matchmean(nmodes::Int, nsamp::Int, E0::Float64, D0::Float64, 
		thup_vec::Vector{Float64}, savemicro::Bool=true)
	# Preliminaries
	nthetas = endof(thup_vec)
	thdn_vec = zeros(Float64,nthetas)
	# Sample H3 and H2 from a microcanonical distribution.
	cputime = @elapsed (H3vec, H2vec, rvar) = microcan(nmodes,nsamp)
	cputime = signif(cputime/60,2)
	println("CPU time for microcanonical sample is ", cputime, " minutes.")
	# For each theta_up, find the corresponding theta_dn by matching the mean.
	meanham_dn(theta_dn::Float64) = meanham(H3vec,H2vec,E0,D0,theta_dn,false)
	for nn = 1:nthetas
		# Determine theta_dn through root finding.
		mean_up = meanham(H3vec,H2vec,E0,D0,thup_vec[nn],true)
		meandiff(theta_dn::Float64) = meanham_dn(theta_dn) - mean_up
		thdn_vec[nn] = find_zero(meandiff, thup_vec[nn], Order1())
		# Sample from the upstream and downstream Gibbs measures
		upsuffix = string("up",nn)
		dnsuffix = string("dn",nn)
		gibbs_sample(rvar,H3vec,H2vec, E0,1.,thup_vec[nn], savemicro,upsuffix)
		gibbs_sample(rvar,H3vec,H2vec, E0,D0,thdn_vec[nn], savemicro,dnsuffix)
	end

end

# TO DO 
# Pass thdn_vec to modify values

#= Main routine to enforce the statistical matching condition. =#
function main(paramsfile::AbstractString="params.txt")
	# Read the parameters from a file.
	params = readvec(paramsfile)
	nmodes, nsampper, npasses = Int(params[1]), Int(params[2]), Int(params[3])
	E0, D0, thmin, thmax, dth = params[3:7]
	thup_vec = collect(thmin:dth:thmax)
	newfolder(datafolder())
	# Take a number of passes, each time with a manageable number of samples. 
	tm0 = time()
	for ii=1:npasses
		matchmean(nmodes, nsampper, E0, D0, thup_vec)
	end
	cputime = signif((time()-tm0)/60, 2)
	nsamptot = sampper*npasses
	println("The total CPU time is ", cputime, " minutes.")
	writebasicdata(params,nsamptot,cputime,thup_vec,thdn_vec)
end

#plt = plot(thup_vec,thdn_vec, xlabel="theta_up",ylabel="theta_dn"); display(plt)


