include("routines.jl")

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
	savemicro = true
	# Set the upstream inverse temperatures to use.
	thup_vec = -.4:0.05:0.1
	thdn_vec = zeros(Float64,0)
	# Sample H3 and H2 from a microcanonical distribution.
	H3vec, H2vec, rvar = microcan(nmodes,nsamptot)
	# Define a function for the downstream mean of the downstream Hamiltonian.
	meanham_dn(theta_dn::Float64) = meanham(H3vec,H2vec,E0,D0,theta_dn,false)
	# For each theta_up, find the corresponding theta_dn by matching the mean.
	for theta_up in thup_vec
		# Compute the upstream mean of the downstream Hamiltonian.
		mean_up = meanham(H3vec,H2vec,E0,D0,theta_up,true)
		# Find the theta_dn via a root find.
		meandiff(theta_dn::Float64) = meanham_dn(theta_dn) - mean_up
		theta_dn = find_zero(meandiff, theta_up, Order1())
		push!(thdn_vec,theta_dn)
		# Sample from the upstream and downstream Gibbs measures
		gibbs_sample(H3vec, H2vec, E0, 1., theta_up, savemicro)
		gibbs_sample(H3vec, H2vec, E0, D0, theta_dn, savemicro)
	end
	# Make plot of theta_dn versus theta_up.
	plt = plot(thup_vec,thdn_vec, xlabel="theta_up",ylabel="theta_dn")
	display(plt)
	return 
end

matchmean(10, 1*10^5, 4., 0.5)

