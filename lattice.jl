# lattice.jl
#
# Defines the spatial lattice; supports 1D and 2D square, open and periodic boundaries
# The main data structure is the bond-spin index array, bond_spin[nBond,2]


function lattice_bond_spins(::Chain{L, Periodic}) where {L}
    bond_spin = lattice_bond_spins(Chain(L; boundary = Fixed))
    return vcat(bond_spin, [L 1])
end

function lattice_bond_spins(::Chain{L, Fixed}) where {L}
    return [i + j for i in 1:(L-1), j in 0:1]
end

# 2D square lattice with fixed/open BCs
function lattice_bond_spins(::Square{Tuple{L,L},Fixed} where {L})
    nSpin = L * L
    nBonds = 2 * nSpin - L - L
    bond_spin = zeros(Int, nBonds, 2)

    # horizontal
    cnt = 1
    @inbounds for i in 1:div(nBonds, 2)
        if mod(i, L - 1) != 0
            bond_spin[i, 1] = cnt
            bond_spin[i, 2] = cnt + 1
            cnt += 1
        else
            bond_spin[i, 1] = cnt
            bond_spin[i, 2] = cnt + 1
            cnt += 2
        end
    end

    # vertical
    cnt = 1
    @inbounds for i in (div(nBonds, 2)+1):nBonds
        bond_spin[i, 1] = cnt
        bond_spin[i, 2] = cnt + L
        cnt += 1
    end

    return bond_spin
end


# 2D square lattice with periodic BCs
function lattice_bond_spins(::Square{Tuple{L,L},Periodic} where {L})
    nSpin = L * L
    nBonds = 2 * nSpin
    bond_spin = zeros(Int, nBonds, 2)  # assign site indices to bonds

    @inbounds for i in 1:nSpin
        # horizontal
        bond1 = 2 * i - 1
        bond_spin[bond1, 1] = i
        bond_spin[bond1, 2] = i + 1

        if mod(i, L) == 0
            bond_spin[bond1, 2] -= L
        end

        # vertical
        bond2 = 2 * i
        bond_spin[bond2, 1] = i
        bond_spin[bond2, 2] = i + L

        if i > (nSpin - L)
            bond_spin[bond2, 2] -= nSpin
        end
    end

    return bond_spin
end
