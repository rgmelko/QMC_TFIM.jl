
abstract type OperatorForm end
struct Diagonal <: OperatorForm end
struct OffDiagonal <: OperatorForm end

abstract type SSEOperator{N, F <: OperatorForm, L <: BoundedLattice} end

struct SiteOperator{F, L} <: SSEOperator{1, F, L}
    i::Int
end
struct IdOperator{L} <: SSEOperator{1, Diagonal, L}
    i::Int
end
struct BondOperator{F, L} <: SSEOperator{2, F, L}
    i::Int
    j::Int
end

function SiteOperator{F}(i::Int, lat::L) where {L <: BoundedLattice, F <: OperatorForm}
    @assert 0 < i <= length(lat)
    return SiteOperator{F, L}(i)
end

function IdOperator(i::Int, lat::L) where L <: BoundedLattice
    @assert 0 < i <= length(lat)
    return IdOperator{L}(i)
end
function BondOperator{F}(i::Int, j::Int, lat::L) where {L <: BoundedLattice, F <: OperatorForm}
    @assert 0 < i <= length(lat)
    @assert 0 < j <= length(lat)
    return BondOperator{F, L}(i, j)
end

operatorform(::Type{<:SSEOperator{N, F}}) where {N, F} = F
operatorform(::T) where {T <: SSEOperator} = operatorform(T)

isdiagonal(T::Type{<:SSEOperator}) = (operatorform(T) <: Diagonal)
isdiagonal(::T) where {T <: SSEOperator} = isdiagonal(T)

issiteoperator(T::Type{<:SSEOperator}) = (T <: SiteOperator)
issiteoperator(::T) where {T <: SSEOperator} = issiteoperator(T)

isbondoperator(T::Type{<:SSEOperator}) = (T <: BondOperator)
isbondoperator(::T) where {T <: SSEOperator} = isbondoperator(T)
