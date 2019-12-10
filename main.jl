# main.jl

include("lattice.jl") #define the spatial lattice

#Projector parameters

M = 10 #length of the projector operator_list is 2M

spin_left = fill(1,nSpin) #left and right trail spin state
spin_right = fill(1,nSpin)

operator_list = fill(0,(2*M,2))
#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j)

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




println(M)
println(spin_left)
println(spin_right)
println(operator_list)

rr = rand()
println(rr)

