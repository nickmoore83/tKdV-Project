## RUN STUFF
include("routines.jl")

# Test case with Lambda = 16, nsweeps = 10, lamfac = 8.
#uniform_sample(16, 2)
#transfun("rand-16z-100.jld", 8)
#output_text("thvars-16z-8-10.jld")

# RUN
#transfun("rand-16z-500.jld", 8)
transfun("rand-16z-1000.jld", 8)
transfun("rand-32z-500.jld", 8)
transfun("rand-32z-1000.jld", 8)
transfun("rand-32z-500.jld", 16)
transfun("rand-32z-1000.jld", 16)