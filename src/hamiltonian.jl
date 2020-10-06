abstract type Hamiltonian{D,N} end

localdim(::Hamiltonian{D}) where {D} = D
dim(::Hamiltonian{D,N}) where {D,N} = N

struct TFIM{N} <: Hamiltonian{2,N}
    op_sampler::OperatorSampler{2, Float64}
    h::Float64
    J::Float64
    P_normalization::Float64
    Ns::Int
    Nb::Int
end

function TFIM(bond_spin, Dim::Int, Ns::Int, Nb::Int, h::Float64, J::Float64)
    ops, p = make_prob_vector(bond_spin, Ns, J, h)
    op_sampler = OperatorSampler(ops, p)
    P_normalization = sum(p)
    return TFIM{Dim}(op_sampler, h, J, P_normalization, Ns, Nb)
end

zero(H::Hamiltonian{2}) = falses(nspins(H))
zero(H::Hamiltonian) = zeros(nspins(H))
one(H::Hamiltonian{2}) = trues(nspins(H))
one(H::Hamiltonian) = ones(nspins(H))

nspins(H::TFIM) = H.Ns
nbonds(H::TFIM) = H.Nb
