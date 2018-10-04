include("routines.jl")

# Simple test for gibbs
function gibbs0test()
	Lambda = 8
	etah = gibbs0(Lambda)
	Ecomp = energy(etah)
	E0 = 1.0
	reldiff = abs((E0-Ecomp)/E0)
	println("\nGibbs0 test")
	println("The desired energy is E0 = ", signif(E0,3))
	println("The energy of etahat is ", signif(Ecomp,3))
	println("The relative difference is ", signif(reldiff,3))
end

#= Test for ham3 and ham3sum.
If eta = 2*a1*cos(x) - 2*b1*sin(x) + 2*a2*cos(2x) - 2*b1*sin(2x)
then eta0 = 0, eta1 = a1 + im*b1, eta2 = a2 + im*b2 
and H3 = a1^2*a2 + 2*a1*b1*b2 - a2*b1^2 =#
function H3test()
	# Parameters
	nmodes = 8
	a1 = 4.3
	b1 = -1.2
	a2 = 2.3
	b2 = 4.5
	# Calculations
	H3exact = a1^2*a2 + 2*a1*b1*b2 - a2*b1^2
	etah = zeros(Complex128,nmodes)
	etah[1] = a1 + im*b1
	etah[2] = a2 + im*b2
	H3comp = ham3(etah)
	H3sum = ham3sum(etah)
	H3err = abs((H3exact-H3comp)/H3exact)
	H3sumerr = abs((H3exact-H3sum)/H3exact)
	# Print statements
	println("\nH3 test")
	println("The exact H3 is ", signif(H3exact,4))
	println("H3 computed by fft to physical space is ", signif(H3comp,4))
	println("H3 computed by direct summation in spectral space is ", signif(H3sum,4))
	println("Relative error for fft: ", signif(H3err,3))
	println("Relative error for direct: ", signif(H3sumerr,3))
	return
end
#= Test for ham2 .=#
function H2test()
	# Parameters
	nmodes = 4
	a1 = 0.1
	b1 = -1.2
	a2 = -0.4
	b2 = 0.8
	# Calculations
	H2exact = a1^2 + b1^2 + 4*a2^2 + 4*b2^2
	etah = zeros(Complex128,nmodes)
	etah[1] = a1 + im*b1
	etah[2] = a2 + im*b2
	H2comp = ham2(etah)
	relerr = abs((H2exact-H2comp)/H2exact)
	# Print statements
	println("\nH2 test")
	println("The exact H2 is ", signif(H2exact,4))
	println("The computed H2 is ", signif(H2comp,4))
	println("The relative error is ", signif(relerr,3))
end
# Test my real fft and its inverse.
function ffttest()
	# Parameters
	nmodes = 32
	npts = npoints(nmodes)
	a0 = 0.0
	a1 = -0.8
	b1 = 1.3
	a2 = -0.2
	b2 = -1.4
	a3 = -0.1
	b3 = 0.8
	# Calculations
	dx = 2*pi/npts
	xx = range(0, dx, npts)
	uu = a0 + 2*a1*cos(xx) - 2*b1*sin(xx) + 2*a2*cos(2*xx) - 2*b2*sin(2*xx)
	uu += 2*a3*cos(3*xx) - 2*b3*sin(3*xx)
	# Calculate uhat.
#	uhat = rfft(uu)/nmodes
	uhat = realfft(uu)
	uback = irealfft(uhat)
	ubackup = ifftup(uhat)
	uhdiffvec = [uhat[1]-(a1+im*b1), uhat[2]-(a2+im*b2), uhat[3]-(a3+im*b3)]
	uhaterr = maxabs(uhdiffvec)/maxabs(uhat)
	ureldiff1 = maxabs(uu-uback)/maxabs(uu)
	#ureldiff2 = maxabs(uu-ubackup)/maxabs(uu)
	# Display vectors.
#	display(round(uhat,6))
#	display(uu)
#	display(uback)
	# Print statements.
	println("\nReal FFT tests")
	println("length of uu ", length(uu))
	println("length of uhat ", length(uhat))
	println("length of uback ", length(uback))
	println("relative error in uhat ", signif(uhaterr,3))
	println("relative difference between original u and returned u ", signif(ureldiff1,3))
	#println("same but with upsampling ", signif(ureldiff2,3))
end

#= Test my guess that the maximum of H3 over E=E0 sphere is given by
the Dirichlet kernel, with H3 = 1/2 E^(3/2) (Lambda^(1/2)-Lambda^(-1/2)) =#
function H3maxtest()
	nsamples = 1*10^5	# 2*10^5 takes about 10-20 secs.
	Lambda = 6
	E0 = 1.
	maxH3guess = 0.5*E0^(3/2)*sqrt(Lambda)*(1.-1./Lambda)
	maxH3 = 0.
	ehmax = 0.
	for nn=1:nsamples
		etah = gibbs0(Lambda)
		H3 = ham3(etah)
		if maxabs(H3) > maxH3
			maxH3 = maxabs(H3)
			ehmax = etah
		end
	end
	# Test the Dirichlet kernel directly.
	uhdir = sqrt(E0/Lambda) * ones(Complex128,Lambda)
	Edir = energy(uhdir)
	H3dir = ham3(uhdir)
	println("\nFirst simply testing the Dirichlet kernel uDir")
	println("The energy of uDir is ", signif(Edir,3))
	println("H3 is ", signif(H3dir,3))
	println("My guess was ", signif(H3dir,3))
	# View the maximum H3 encountered.
	println("\nTest for the max of H3")
	println("The max of H3 encountered is ", signif(maxH3,3))
	println("My guess for max of H3 is ", signif(maxH3guess,3))
	println("\nThe coefficients corresponding to the max are:")
	display(ehmax)
	println("The energy of each mode is:")
	display(abs(ehmax).^2)
	# Plot the corresponding eta
	etamax = ifftup(ehmax)
	npts = endof(etamax)
	dx = 2*pi/npts
	xx = range(0,dx,npts)
	p1 = plot([xx; 2*pi], [etamax; etamax[1]])
	display(p1)
end

#= Compare CPU time for sampling all at once versus one at a time. =#
function sampletest()
	# Take one sample.
	function sampleone(nmodes::Int)
		rvar = randn(nmodes,2)
		uhat = rvar[:,1] + im*rvar[:,2]
		return uhat
	end
	# Take all the samples at once.
	function sampleall(nmodes::Int,nsamples::Int)
		uhat = zeros(Complex128,nmodes)
		H3 = zeros(Float64,nsamples)
		H2 = zeros(Float64,nsamples)
		rvar = randn(nmodes,2,nsamples)
		for nn=1:nsamples
			uhat = rvar[:,1,nn] + im*rvar[:,2,nn]
			H3[nn] = ham3(uhat)
			H2[nn] = ham2(uhat)
		end
		return H3,H2
	end
	# Test the two methods.
	nmodes = 10
	nsamples = 2*10^5
	# Sample one at a time.
	println("Sampling one at a time.")
	cput0 = time()
	H3one = zeros(Float64,nsamples)
	H2one = zeros(Float64,nsamples)
	for nn=1:nsamples
		uhat = sampleone(nmodes)
		H3one[nn] = ham3(uhat)
		H2one[nn] = ham2(uhat)
	end
	cputimeone = time()-cput0
	# Sample all at once.
	println("Sampling all at once.")
	cput0 = time()
	H3all,H2all = sampleall(nmodes,nsamples)
	cputimeall = time()-cput0
	# Print results.
	println("\nThe size of H3all is ", size(H3all))
	println("The size of H3one is ", size(H3one))
	println(); display(H3all)
	println(); display(H3one)
	println("\nCPU time to sample all at once is ", signif(cputimeall,3))
	println("CPU time to sample one at a time is ", signif(cputimeone,3))
end

#= Test breaking out of a loop. =#
function testbreak()
	ii = 0
	for ii=1:10
		println("In loop, ii is ", ii)
		if ii >=5
			break
		end
	end
	println("Out of loop, ii is ", ii)
end

#= Test parallel coding =#
function parmain()
	nmodes = 16
	nsamples = 10^8
	# Get all the random samples and allocate space.
	#rvar = randn(nmodes,2,nsamples)
	#H3vec, H2vec = [zeros(Float64,nsamples) for nn=1:2]
	rvar = SharedArray{Float64}( randn(nmodes,2,nsamples) )
	H3vec, H2vec = [SharedVector{Float64}(nsamples) for nn=1:2]
	@parallel for nn=1:nsamples
		uhat = getuhat(rvar,nn)
		H3vec[nn] = ham3(uhat)
		H2vec[nn] = ham2(uhat)
		if mod(nn, 10^4) == 0
			println("Microcanonical sampling is ", signif(100*nn/nsamples,3), "% completed.")
		end
	end
	println("\nOut of for loop.")
	return RandSet(H3vec,H2vec,rvar)
end

rset = parmain()


# TO DO 
# Also test scaling prediciton for energy of u_x.

#gibbs0test()
#H3test()
#H2test()
#ffttest()

#H3maxtest()
#sampletest()




