# lattice.jl
#
# Defines the spatial lattice; supports 1D and 2D square, open and periodic boundaries
# The main data structure is the bond-spin index array, bond_spin[nBond,2]

using Random
Random.seed!(1234)

Dim = 1
nX = 4
PBC = true

if Dim == 2
	nY = nX  #Works for square lattices only right now
end

#one-dimensional lattice
if Dim == 1
    nSpin = nX
    if PBC == true
       nBond = nSpin
    else
       nBond = nSpin-1
    end

    bond_spin = zeros(Int,nBond,2) #assign site indices to bonds
    for i = 1:nBond
        bond_spin[i,1] = i
        bond_spin[i,2] = i+1
    end
    if PBC == true
        bond_spin[nBond,2] = 1
    end

    #println(bond_spin)

#two-dimensional square lattice OBC
elseif Dim == 2 && PBC == false
    nSpin = nX*nY
    nBond = 2*nSpin-nX-nY
    bond_spin = zeros(Int,nBond,2) #assign site indices to bonds

    cnt=1
    for i = 1:div(nBond,2) #horizontal
        #println(i," ",cnt)
        if mod(i,nX-1) != 0 
            bond_spin[i,1] = cnt 
            bond_spin[i,2] = cnt + 1
            global cnt += 1
        else
            bond_spin[i,1] = cnt 
            bond_spin[i,2] = cnt + 1
            global cnt += 2
        end
    end
    global cnt = 1
    for i = (div(nBond,2)+1):nBond #vertical
        bond_spin[i,1] = cnt 
        bond_spin[i,2] = cnt + nX
        global cnt += 1
    end

#two-dimensional square lattice PBC
elseif Dim == 2 && PBC == true
    nSpin = nX*nY
    nBond = 2*nSpin
    bond_spin = zeros(Int,nBond,2) #assign site indices to bonds

    for i = 1:nSpin
        bond1 = 2*i - 1 #horizontal
        bond_spin[bond1,1] = i
        bond_spin[bond1,2] = i+1
        if mod(i,nX) == 0
            bond_spin[bond1,2] -= nX
        end

        bond2 = 2*i     #vertical
        bond_spin[bond2,1] = i
        bond_spin[bond2,2] = i+nX
        if i > (nSpin - nX)
            bond_spin[bond2,2] -= nSpin
        end
    end

else 
    println("Dimension error")
end

# println("Dimension ",Dim)
# println("PBC ",PBC)
# println("Number of bonds ", nBond)
# println(bond_spin)
# 
# for i = 1:nBond
#     println(i," ",bond_spin[i,1]," ",bond_spin[i,2])
# end
