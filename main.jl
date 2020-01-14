# main.jl
#
# A projector QMC program for the TFIM

using Random
Random.seed!(1234)

include("lattice.jl") #define the spatial lattice
include("updates.jl") #functions for the Monte Carlo updates
include("measurements.jl") #functions for the Monte Carlo updates


Dim = 1
nX = 6
PBC = false

lattice_shape = tuple(repeat([nX], Dim)...)
BC = PBC ? Periodic : Fixed

lattice = HyperRect(lattice_shape; boundary = BC)
bond_spin = lattice_bond_spins(lattice)

nSpin = prod(size(lattice))
nBond = size(bond_spin, 1)

#Projector parameters
M = 300 #length of the projector operator_list is 2M
h_x = 0.5
J_ = 1.0
MCS = 10000 #the number of Monte Carlo steps

#*******  Globals
spin_left = fill(1, nSpin) #left and right trail spin state
spin_right = fill(1, nSpin)
LinkList = []
LegType = []
Associates = []
operator_list = fill(0, (2 * M, 2))
#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j)
#*******  Globals

for i in 1:2000  #Equilibration
    DiagonalUpdate()
    LinkedList()
    ClusterUpdate()
end

M2 = 0
for i in 1:MCS #Monte Carlo Production Steps
    DiagonalUpdate()
    LinkedList()
    global M2 += Measure()
    ClusterUpdate()
end

println(M2 / (nSpin * nSpin * MCS))
