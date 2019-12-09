# main.jl

include("lattice.jl") #define the spatial lattice

#include projector
struct Basis 
    M::Int  #length of the projector list is 2M
	spin_left::Vector{Int}  #left and right trail spin state
	spin_right::Vector{Int}
end#struct

for i=1:1000
    rr = rand()
    println(rr)
end

rspin = rand([0,1],nSpin) #Random spin configuration
B = Basis(3,rspin,rspin)

println(B.M)
println(B.spin_left)
println(B.spin_right)
