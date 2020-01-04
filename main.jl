# main.jl
#
# A projector QMC program for the TFIM

include("lattice.jl") #define the spatial lattice
include("updates.jl") #functions for the Monte Carlo updates
include("measurements.jl") #functions for the Monte Carlo updates

#Projector parameters
M = 300 #length of the projector operator_list is 2M
h_x = 0.5
J_ = 1.0
MCS = 10000 #the number of Monte Carlo steps

#*******  Globals
spin_left = fill(1,nSpin) #left and right trail spin state
spin_right = fill(1,nSpin)
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


for i = 1:2000  #Equilibration
    DiagonalUpdate()
    LinkedList()
    ClusterUpdate()
end

M2 = 0
for i = 1:MCS #Monte Carlo Production Steps
    DiagonalUpdate()
    LinkedList()
	global M2 += Measure()
    ClusterUpdate()
end

println(M2/(nSpin*nSpin*MCS))
