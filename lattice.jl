# lattice.jl
#
# Defines the spatial lattice; supports 1D and 2D square, open and periodic boundaries
# The main data structure is the bond-spin index array, bond_spin[nBond,2]


@inline lattice_bond_spins(C::Chain{L}) where {L} = collect(edges(C))

# 2D square lattice with fixed/open BCs
function lattice_bond_spins(::Square{Tuple{L,L},Fixed}) where {L}
    nSpin = L * L
    nBonds = 2 * nSpin - L - L
    bond_spin = [(0, 0) for _ in 1:nBonds]

    # horizontal
    cnt = 1
    @inbounds for i in 1:div(nBonds, 2)
        if mod(i, L - 1) != 0
            bond_spin[i] = (cnt, cnt + 1)
            cnt += 1
        else
            bond_spin[i] = (cnt, cnt + 1)
            cnt += 2
        end
    end

    # vertical
    cnt = 1
    @inbounds for i in (div(nBonds, 2)+1):nBonds
        bond_spin[i] = (cnt, cnt + L)
        cnt += 1
    end

    return bond_spin
end


# 2D square lattice with periodic BCs
function lattice_bond_spins(::Square{Tuple{L,L},Periodic}) where {L}
    nSpin = L * L
    nBonds = 2 * nSpin
    bond_spin = [(0, 0) for _ in 1:nBonds]  # assign site indices to bonds

    @inbounds for i in 1:nSpin
        # horizontal
        bond1 = 2 * i - 1

        if mod(i, L) != 0
            bond_spin[bond1] = (i, i + 1)
        else
            bond_spin[bond1] = (i, i + 1 - L)
        end

        # vertical
        bond2 = 2 * i

        if i <= (nSpin - L)
            bond_spin[bond2] = (i, i + L)
        else
            bond_spin[bond2] = (i, i + L - nSpin)
        end
    end

    return bond_spin
end
