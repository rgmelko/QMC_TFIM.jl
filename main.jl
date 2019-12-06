#using Random
#rng = MersenneTwister(1111)

Dim = 1
nX = 4
nY = 4
PBC = true

#Spin = Array{Int,Dim}[]
#S = rand(rng)

Spin = rand([0,1],nX,nY)
Spin = 2*Spin - 1

println(Dim)
println(PBC)
println(Spin)
