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
