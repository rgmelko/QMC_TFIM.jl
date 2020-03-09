import Base: zero


include("lattice.jl") # define the spatial lattice


function init_op_list(length)
    operator_list::Vector{NTuple{2,Int}} = [(0, 0) for _ in 1:length]
    return operator_list
end

function grow_op_list!(operator_list::Vector{NTuple{2,Int}}, factor::Real)
    len = round(Int, (factor - 1) * length(operator_list))
    tail = init_op_list(len)
    append!(operator_list, tail)
end

abstract type Hamiltonian{D,N,L<:BoundedLattice{N}} end

localdim(::Hamiltonian{D}) where {D} = D
dim(::Hamiltonian{D,N}) where {D,N} = N

struct TFIM{N,L} <: Hamiltonian{2,N,L}
    lattice::L
    bond_spin::Vector{NTuple{2, Int}}
    h::Float64
    J::Float64
    P_h::Float64
    P_J::Float64
    P_normalization::Float64
    Ns::Int
    Nb::Int
end

function TFIM(lattice::L, h::Float64, J::Float64) where {L<:BoundedLattice{N}} where {N}
    bond_spin = lattice_bond_spins(lattice)

    Ns, Nb = length(lattice), length(bond_spin)

    # define the Metropolis probability as a constant
    # https://pitp.phas.ubc.ca/confs/sherbrooke2012/archives/Melko_SSEQMC.pdf
    # equation 1.43
    P_h = h * Ns
    P_J = 2 * J * Nb

    P_normalization = P_h + P_J
    P_h /= P_normalization
    P_J /= P_normalization

    return TFIM{N,L}(lattice, bond_spin, h, J, P_h, P_J, P_normalization, Ns, Nb)
end

zero(H::Hamiltonian{2}) = falses(size(H.lattice)...)
zero(H::Hamiltonian) = zeros(size(H.lattice)...)
one(H::Hamiltonian{2}) = trues(size(H.lattice)...)
one(H::Hamiltonian) = ones(size(H.lattice)...)

nspins(H::TFIM) = H.Ns
nbonds(H::TFIM) = H.Nb

abstract type AbstractQMCState{D,N,H<:Hamiltonian{D,N}} end

struct BinaryQMCState{N,H} <: AbstractQMCState{2,N,H}
    left_config::BitArray{N}
    right_config::BitArray{N}
    operator_list::Vector{NTuple{2,Int}}
end

function BinaryQMCState(H::Hamiltonian{2,N}, M::Int) where {N}
    BinaryQMCState{N,typeof(H)}(zero(H), zero(H), init_op_list(2*M))
end

struct PottsQMCState{D,N,H} <: AbstractQMCState{D,N,H}
    left_config::Array{Int,N}
    right_config::Array{Int,N}
    operator_list::Vector{NTuple{2,Int}}
end

function PottsQMCState(H::Hamiltonian{D,N}, M::Int) where {D,N}
    PottsQMCState{N,typeof(H)}(zero(H), zero(H), init_op_list(2*M))
end
