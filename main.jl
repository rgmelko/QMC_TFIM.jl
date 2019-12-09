using Random
Random.seed!(1234)

Dim = 2
nX = 4
nY = 4
PBC = true

#one-dimensional lattice
if Dim == 1
    nSpin = nX
    Spin = rand([0,1],nX) #Random spin configuration
    if PBC == true
       nBond = nSpin
    else
       nBond = nSpin-1
    end

    bond_spin = zeros(nBond,2) #assign site indices to bonds
    for i = 1:nBond
        bond_spin[i,1] = i
        bond_spin[i,2] = i+1
    end
    if PBC == true
        bond_spin[nBond,2] = 1
    end

    println(bond_spin)

#two-dimensional square lattice OBC
elseif Dim == 2 && PBC == false
    nSpin = nX*nY
    Spin = rand([0,1],nX,nY) #Random spin configuration
    if PBC == true
       nBond = 2*nSpin
    else
       nBond = 2*nSpin-nX-nY
    end

    bond_spin = zeros(nBond,2) #assign site indices to bonds
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
    println("PBC 2D ")

#    if PBC == true
#        bond_spin[nBond,2] = 1
#    end

else 
    println("Dimension error")
end

println("Dimension ",Dim)
println("PBC ",PBC)
println("Spin config ",Spin)
println("Number of bonds ", nBond)
println(bond_spin)

for i = 1:nBond
    println(i," ",bond_spin[i,1]," ",bond_spin[i,2])
end
