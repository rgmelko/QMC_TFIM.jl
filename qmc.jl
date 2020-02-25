import Base: zero


include("lattice.jl") # define the spatial lattice


function init_op_list(lattice, M)
    operator_list::Vector{NTuple{2,Int}} = [(0, 0) for _ in 1:2M]
    return operator_list
end

abstract type Hamiltonian{D,N,L<:BoundedLattice{N}} end

localdim(::Hamiltonian{D}) where {D} = D
dim(::Hamiltonian{D,N}) where {D,N} = N

struct TFIM{N,L} <: Hamiltonian{2,N,L}
    lattice::L
    bond_spin::Vector{NTuple{2, Int}}
    h::Float64
    J::Float64
end

function TFIM(lattice::L, h::Float64, J::Float64) where {L<:BoundedLattice{N}} where {N}
    bond_spin = lattice_bond_spins(lattice)
    return TFIM{N,L}(lattice, bond_spin, h, J)
end

zero(H::Hamiltonian{2}) = falses(size(H.lattice)...)
zero(H::Hamiltonian) = zeros(size(H.lattice)...)
one(H::Hamiltonian{2}) = trues(size(H.lattice)...)
one(H::Hamiltonian) = ones(size(H.lattice)...)

nspins(H::TFIM) = length(H.lattice)
nbonds(H::TFIM) = size(H.bond_spin, 1)

abstract type AbstractQMCState{D,N,H<:Hamiltonian{D,N}} end

struct BinaryQMCState{N,H} <: AbstractQMCState{2,N,H}
    left_config::BitArray{N}
    right_config::BitArray{N}
    operator_list::Vector{NTuple{2,Int}}
end

function BinaryQMCState(H::Hamiltonian{2,N}, M::Int) where {N}
    BinaryQMCState{N,typeof(H)}(zero(H), zero(H), init_op_list(H.lattice, M))
end

struct PottsQMCState{D,N,H} <: AbstractQMCState{D,N,H}
    left_config::Array{Int,N}
    right_config::Array{Int,N}
    operator_list::Vector{NTuple{2,Int}}
end

function PottsQMCState(H::Hamiltonian{D,N}, M::Int) where {D,N}
    PottsQMCState{N,typeof(H)}(zero(H), zero(H), init_op_list(H.lattice, M))
end
