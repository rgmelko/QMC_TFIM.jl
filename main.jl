# main.jl
#
# A projector QMC program for the TFIM

include("lattice.jl") #define the spatial lattice
include("updates.jl") #functions for the Monte Carlo updates

#Projector parameters
M = 4 #length of the projector operator_list is 2M
h_x = 1.0
J_ = 1.0

#define the Metropolis probability as a constant
#https://pitp.phas.ubc.ca/confs/sherbrooke2012/archives/Melko_SSEQMC.pdf
#equation 1.43
const P_h = h_x*nSpin/(h_x*nSpin +2.0*J_*nBond) #J=1.0 tested only

#*******  Globals
spin_left = fill(1,nSpin) #left and right trail spin state
spin_right = fill(1,nSpin)
println(spin_left)

LinkList =[]
LegType =[]
Associates =[]

operator_list = fill(0,(2*M,2))
#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j)

#*******  Globals

#initialize the operator list
for i = 1:2*M
    coin = rand(Bool)
    if coin == true #insert a diagonal site operator
        site = rand(1:nSpin)
        operator_list[i,1] = -1
        operator_list[i,2] = site
    else
        bond = rand(1:nBond)
        operator_list[i,1] = bond_spin[bond,1]
        operator_list[i,2] = bond_spin[bond,2]
    end
end

DiagonalUpdate()
LinkedList()
ClusterUpdate()
#println(operator_list)
#println(spin_left)
#println(spin_right)

DiagonalUpdate()
LinkedList()
ClusterUpdate()
println(operator_list)
println(spin_left)
println(spin_right)

#propagate the spin state through the operator list
#spin_prop = copy(spin_left)
#for i = 1:2*M
#    #spin_prop[1] = 0
#end

#println(M)
#println(spin_left)
#println(spin_right)
#println(spin_prop)


