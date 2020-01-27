
abstract type OperatorForm end
struct Diagonal <: OperatorForm end
struct OffDiagonal <: OperatorForm end

# N is the number of sites the operator acts on
abstract type SSEOperator{N, F <: OperatorForm} end

struct SiteOperator{F} <: SSEOperator{1, F}
    i::Int
end
struct IdOperator <: SSEOperator{1, Diagonal}
    i::Int
end
struct BondOperator{F} <: SSEOperator{2, F}
    i::Int
    j::Int
end

function SiteOperator{F}(i::Int, lat::L) where {L <: BoundedLattice, F <: OperatorForm}
    @assert 0 < i <= length(lat)
    return SiteOperator{F}(i)
end

function IdOperator(i::Int, lat::L) where L <: BoundedLattice
    @assert 0 < i <= length(lat)
    return IdOperator(i)
end

function BondOperator{F}(i::Int, j::Int, lat::L) where {L <: BoundedLattice, F <: OperatorForm}
    @assert 0 < i <= length(lat)
    @assert 0 < j <= length(lat)
    return BondOperator{F}(i, j)
end

operatorform(::Type{<:SSEOperator{N, F}}) where {N, F} = F
operatorform(::T) where {T <: SSEOperator} = operatorform(T)

isdiagonal(T::Type{<:SSEOperator}) = (operatorform(T) <: Diagonal)
isdiagonal(::T) where {T <: SSEOperator} = isdiagonal(T)

issiteoperator(T::Type{<:SSEOperator}) = (T <: SiteOperator)
issiteoperator(::T) where {T <: SSEOperator} = issiteoperator(T)

isbondoperator(T::Type{<:SSEOperator}) = (T <: BondOperator)
isbondoperator(::T) where {T <: SSEOperator} = isbondoperator(T)
