using Random
Random.seed!(1234)

Dim = 1
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

#two-dimensional lattice
elseif Dim == 2
    nSpin = nX*nY
    Spin = rand([0,1],nX,nY) #Random spin configuration
    if PBC == true
       nBond = 2*nSpin
    else
       nBond = nSpin-nX-nY
    end

else 
    println("Dimension error")
end

println("Dimension ",Dim)
println("PBC ",PBC)
println("Spin config ",Spin)
println("Number of bonds ", nBond)
