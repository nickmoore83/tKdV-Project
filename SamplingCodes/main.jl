include("routines.jl")

#= Write the basic data to a file. =#
function writebasicdata(nmodes::Int, nsamptot::Int, E0::Float64, D0::Float64,
		cputime::Float64, thup_vec::Vector{Float64}, thdn_vec::Vector{Float64})
	foldername = datafolder()
	file = string(foldername,"basic.txt")
	label1 = "# Basic parameters: nmodes, nsamptot, E0, D0, CPU time (mins)"
	label2 = "# Inverse temperature data: number of thetas, theta_ups and theta_dns"
	data = [label1; nmodes; nsamptot; E0; D0; cputime; 
		label2; endof(thup_vec); thup_vec; thdn_vec]
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

#= Enforce the statistical matching condition. =#
function matchmean(nmodes::Int, nsamptot::Int, E0::Float64, D0::Float64)
	# Preliminaries
	thup_vec = collect(-.4:0.05:0.1)
	#thup_vec = [-0.1]
	savemicro = true
	newfolder(datafolder())
	nthetas = endof(thup_vec)
	thdn_vec = zeros(Float64,nthetas)
	# Sample H3 and H2 from a microcanonical distribution.
	cputime = @elapsed (H3vec, H2vec, rvar) = microcan(nmodes,nsamptot)
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
	writebasicdata(nmodes,nsamptot,E0,D0,cputime,thup_vec,thdn_vec)
	#plt = plot(thup_vec,thdn_vec, xlabel="theta_up",ylabel="theta_dn"); display(plt)
	return 
end

# Quick testing
#matchmean(10, 1*10^5, 4., 0.5)

matchmean(20, 1*10^7, 4., 0.5)

#matchmean(20, 4*10^7, 4., 0.5)

