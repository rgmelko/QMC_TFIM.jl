# main.jl

include("lattice.jl") #define the spatial lattice

rng = MersenneTwister(1133);

for i=1:1000
    rr = rand(rng,Float64)
    println(rr)
end
