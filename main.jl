# main.jl
#
# A projector QMC program for the TFIM

using Random
Random.seed!(1234)

using DelimitedFiles

include("lattice.jl") #define the spatial lattice
include("updates.jl") #functions for the Monte Carlo updates
include("measurements.jl") #functions for the Monte Carlo updates


Dim = 1
nX = 10
PBC = false

lattice_shape = tuple(repeat([nX], Dim)...)
BC = PBC ? Periodic : Fixed
BC_name = PBC ? "PBC" : "OBC"

lattice = HyperRect(lattice_shape; boundary = BC)
bond_spin = lattice_bond_spins(lattice)

nSpin = prod(size(lattice))
nBond = size(bond_spin, 1)

#Projector parameters
M = 300 #length of the projector operator_list is 2M
h = 1.0
J_ = 1.0
MCS = 10000 #the number of Monte Carlo steps

root = "./data/$(Dim)D/$(nX)/$(BC_name)/h$(h)/"
mkpath(root)
output_file = "$(root)samples.txt"

#*******  Globals
spin_left = falses(nSpin) #left and right trail spin state
spin_right = falses(nSpin)
LinkList = []
LegType = []
Associates = []
operator_list = zeros(Int, 2 * M, 2)
#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j + 1)
#*******  Globals

for i in 1:2000  #Equilibration
    diagonal_update!(operator_list, h, J_, nSpin, nBond, spin_left, spin_right)
    LinkedList()
    ClusterUpdate()
end

measurements = falses(MCS, nSpin)

mag, mag2 = 0., 0.
energy = 0.0
for i in 1:MCS #Monte Carlo Production Steps
    diagonal_update!(operator_list, h, J_, nSpin, nBond, spin_left, spin_right)
    LinkedList()
    spin_prop = sample(spin_left, operator_list)

    measurements[i, :] = spin_prop

    global energy += energy_abs_zero(h, J_, spin_prop, operator_list)
    global mag += magnetization(spin_prop)
    global mag2 += magnetization_sqr(spin_prop)
    ClusterUpdate()
end


open(output_file, "a") do io
    writedlm(io, measurements, " ")
end

println(mag / (nSpin * MCS))
println(mag2 / (nSpin * nSpin * MCS))

println(energy / (nSpin * MCS))