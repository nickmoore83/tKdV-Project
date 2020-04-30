## RUN STUFF
include("routines.jl")

# Test case with Lambda = 16, nsweeps = 2, lamfac = 8.
# Could also use @elapsed() but cputime is printed and stored in file anyway.
#uniform_sample(16, 2)
#transfun("rand-16z-2.jld", 8)
#output_text("thvars-16z-8-2.jld")


# RUNS
uniform_sample(16, 2)
#transfun("rand-12z-1000.jld", 8)
