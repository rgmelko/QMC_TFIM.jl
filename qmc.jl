
include("lattice.jl") # define the spatial lattice


function init_op_list(M)
    operator_list::Vector{SSEOperator} = [IdOperator(1, lattice) for _ in 1:2M]
    return operator_list
end

operator_list = init_op_list(M)

struct QMCState{N, L <: BoundedLattice{N}}
    lattice::L
    bond_spin::Array{Int, 2}
    nspins::Int
    nbonds::Int
    left_config::BitArray{N}
    right_config::BitArray{N}
end

function QMCState(lattice::BoundedLattice)
    bond_spin = lattice_bond_spins(lattice)
    left_config = falses(size(lattice)...)
    right_config = falses(size(lattice)...)
    nspins = length(lattice)
    nbonds = size(bond_spin, 1)

    N = ndims(lattice)
    L = typeof(lattice)

    return QMCState{N, L}(lattice, bond_spin, nspins, nbonds, left_config, right_config)
end
