include("routines.jl")

# The relative error between two numbers or vectors.
function rel_err(vnum, vex, sigdig::Int=3)
	relerr = maximum(abs, vnum - vex) / maximum(abs, vex)
	return sig(relerr,sigdig)
end
# Transform from physical to spectral space; only used in benchmark of H3.
function realfft(uu::Vector{Float64}, symint::Bool=true)
	uhat = rfft(uu)/length(uu)
	@assert(abs(uhat[1])/maximum(abs,uhat) < 1e-6)
	uhat[end] *= 0.5
	uhat = uhat[2:end]
	if(symint) uhat[1:2:end] .*= -1 end
	return uhat
end

#---------- REAL FFT TEST ----------#
# Test the real fft, both forward and backwards. All tests PASSED!
# This code has a lot of nuance and hard-fought wisdom, so keep it around.
function test_realfft(nmodes::Int)
	#---------- PARAMETERS ----------#
	# The interval is either symmetric [-pi,pi] or normal [0,2*pi]
	# The signal is defined by either a trig or complex-exp series. 
	# I believe setting npts = 2*nmodes is the *proper* way to use FFT.
	symint = true
	trig = false
	npts = 2*nmodes
	#---------- THE COMPUTATION ----------#
	# Define the interval as either symmetric or normal.
	dx = 2*pi/npts
	xx = range(0, step = dx, length = npts)
	xup = range(0, step = 0.5*dx, length = 2*npts)
	if symint
		xx = xx .- pi
		xup = xup .- pi
	end
	# Define the phsyical signal with either a trig or complex-exp series
	if trig
		# Random trig coefficients. Note the last bn must be zero!
		an, bn = randn(nmodes), randn(nmodes)
		bn[end] = 0
		uhex = 0.5*(an - im*bn)
		uuex = sum( an[nn]*cos.(nn*xx) + bn[nn]*sin.(nn*xx) for nn=1:nmodes )
		upex = sum( an[nn]*cos.(nn*xup) + bn[nn]*sin.(nn*xup) for nn=1:nmodes )
	else
		# Random complex coefficients. Note the last mode must be real!
		uhex = randn(ComplexF64, nmodes)
		uhex[end] = randn(Float64)
		uuex = 2*real( sum( uhex[kk] * exp.(im*kk*xx) for kk = 1:nmodes) )
		upex = 2*real( sum( uhex[kk] * exp.(im*kk*xup) for kk = 1:nmodes) )
	end
	uhnum = realfft(uuex, symint)
	#uunum = irealfft(uhex, symint)
	plan = makeplan(nmodes)
	upnum = irealfft(uhex, plan, symint)
	# Function to print the coefficients.
	function printuh(idx)
		println("idx = ", idx)
		println("input uh ", sig(uhex[idx],3))
		println("output uh ", sig(uhnum[idx],3))
	end
	# Print the error.
	#[printuh(kk) for kk=1:nmodes]
	# Forward FFT
	println("\nTesting forward real FFT")
	println("Relative error in uhat without last mode: ", rel_err(uhnum[1:end-1],uhex[1:end-1]))
	println("Relative error in uhat with last mode: ", rel_err(uhnum,uhex))
	# Inverse FFT
	#println("\nTesting inverse real FFT")
	#println("Relative error in uu: ", rel_err(uunum,uuex))
	# Upsampled inverse FFT
	println("\nTesting upsampled inverse real FFT")
	println("Relative error in uup: ", rel_err(upnum,upex))
end
#test_realfft(16)
#---------------------------------------#


#---------- HAMILTONIAN ACCURACY TESTS ----------#
#= Compute the exact H2 and H3 for a specific cos^2 function. =#
function ham_ex1(nmodes::Int, aa::Float64, x0::Float64)
	# Define the function.
	npts = 2*nmodes
	xx = -pi: 2*pi/npts : pi-2*pi/npts
	bb = sqrt(8/(pi*(1+16*aa^2)))
	uu = bb*( (aa .+ cos.(xx.-x0)).^2 .- (0.5+aa^2) )
	uhat = realfft(uu)
	# Calculate the exact Hamiltonian values.
	h2ex = 1 + 3/(1+16*aa^2)
	h3ex = 8*sqrt(2)*pi*aa^2 * (16*pi*aa^2 + pi)^(-3/2)
	#= Can also compute the exact uhat in case useful. =#
	if true
		uhex = zeros(ComplexF64, nmodes)
		uhex[1] = aa*bb*exp(-im*x0) 
		uhex[2] = bb/4*exp(-2*im*x0)
		println("Within ham_ex1")
		println("Rel error in uhat: ", rel_err(uhat, uhex))
	end
	return uhat, h2ex, h3ex
end
#= Compute the exact H2 and H3 for a general function with nmodes=2. =#
function ham_ex2(nmodes::Int, alpha::Float64, th1::Float64, th2::Float64)
	R1 = 1/sqrt(2*pi) * cos(alpha)
	R2 = 1/sqrt(2*pi) * sin(alpha)
	# The values of uhat.
	uhat = zeros(ComplexF64, nmodes)
	uhat[1] = R1*exp(im*th1)
	uhat[2] = R2*exp(im*th2)
	# Exact values of H2 and H3.
	h3ex = 2*pi*R2*R1^2*cos(2*th1-th2)
	h2ex = 2*pi*(R1^2 + 4*R2^2)
	return uhat, h2ex, h3ex
end
#= Test the computation of H2 and H3 using both exact solutions above. 
All tests PASSED! Keep the code around. =#
function test_ham(nmodes::Int; method::Int=3)	
	# random coefficients.
	aa = randn()
	x0, th1, th2 = 2*pi*rand(3)
	alpha = 0.5*pi*rand()
	# Compute the exact H2 and H3 using method_1 or method_2
	if method == 1
		uhat, h2ex,h3ex = ham_ex1(nmodes, aa, x0)
	elseif method == 2
		uhat, h2ex, h3ex = ham_ex2(nmodes, alpha, th1, th2)
	# Or simply compute the numerical numerical valures of H3 and compare.
	elseif method == 3
		uhat = randn(ComplexF64, nmodes)
		#uhat[end] = randn(Float64)	
		# The above line is not needed to get agreement between h3fft and h3rec!!!
	else return end
	# Compute the numerical H2 and H3 and the relative errors.
	h2num = ham2(uhat)
	h3rec = ham3rec(uhat)
	plan = makeplan(nmodes)
	h3fft = ham3fft(uhat,plan)
	if method in (1,2)
		# H2 error
		println("\nH2 test")
		println("h2num = ", sig(h2num,3))
		println("h2ex = ", sig(h2ex,3))
		println("Rel error = ", rel_err(h2num, h2ex))
		# H3 error
		println("\nH3 test")
		println("h3fft = ", sig(h3fft,3))
		println("h3rec = ", sig(h3rec,3))
		println("h3ex = ", sig(h3ex,3))
		println("Rel error fft = ", rel_err(h3fft, h3ex))
		println("Rel error rec = ", rel_err(h3rec, h3ex))
	elseif method == 3
		# H3 difference
		println("\nH3 test")
		println("h3fft = ", sig(h3fft,3))
		println("h3rec = ", sig(h3rec,3))
		println("Rel diff = ", rel_err(h3fft, h3rec))
	end
end
#test_ham(16; method=3)


#---------- H3 SPEED TESTS ----------#
# Test the speed of h3fft versus h3rec.
function ham3speed(nmodes::Int, ncalls::Int)
	uhat = randn(ComplexF64, nmodes, ncalls)
	# Make a plan for FFT.
	plan = makeplan(nmodes)
	# Call each routine many times and measure the time.
	tm_fft = @elapsed([h3fft = ham3fft(uhat[:,nn], plan) for nn=1:ncalls])
	tm_rec = @elapsed([h3rec = ham3rec(uhat[:,nn]) for nn=1:ncalls])
	tm_h2 = @elapsed([h2 = ham2(uhat[:,nn]) for nn=1:ncalls])
	#Print computational times
	println("\nTesting speed for nmodes = ", nmodes)
	println("Time for h2: ", sig(tm_h2,3), "secs")
	println("Time for fft: ", sig(tm_fft,3), "secs")
	println("Time for rec: ", sig(tm_rec,3), "secs")
	println("FFT is faster by factor of ", sig(tm_rec/tm_fft,2))	
	# Print values of h3 just to make sure things are being computed correctly.
	#println("\nh3fft = ", sig(h3fft, 3), "; h3rec = ", sig(h3rec,3))
	#println("Rel err: ", rel_err(h3fft, h3rec))
	return tm_fft, tm_rec, tm_h2
end
#ham3speed(64, 10^5)

# Test the speed of h3fft versus h3rec for several nmodes.
function test_ham_speed()
	ncalls = 2*10^3
	#nmodes = 8:1:128
	nmodes = [8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256] #382, 512
	tfft, trec, th2 = [zeros(Float64, length(nmodes)) for ii=1:3]
	for (idx,nm) in enumerate(nmodes)
		tfft[idx], trec[idx], th2[idx] = ham3speed(nm, ncalls)
	end
	ratio = trec./tfft
	# Plot the CPU times for each.
	p1 = plot(nmodes, trec, label="rec", markershape=:circle)
	plot!(p1, nmodes, tfft, label="fft", markershape=:circle)
	plot!(p1, nmodes, th2, label="h2", markershape=:circle)
	plot!(xlabel="Lambda", ylabel="CPU time (sec)")
	plot!(xscale=:lin, yscale=:lin, legend=:topleft)
	# Plot the ratio of CPU times.
	p2 = plot(nmodes, ratio, label="ratio", markershape=:circle)
	plot!(p2, xlabel="Lambda", ylabel="ratio fft/rec", xscale=:lin, yscale=:log)
	# Display
	display(p1)
end
test_ham_speed()


