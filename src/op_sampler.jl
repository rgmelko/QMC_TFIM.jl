# (-3, i) = long field
# (-2, i) = transverse field off-diag
# (-1, i) = transverse field diag
# (0, 0) = id
# (i, j) = diag bond op

function make_prob_vector(J::AbstractMatrix, hx::AbstractVector)
    ops = Vector{NTuple{2, Int}}(undef, 0)
    p = Vector{Float64}(undef, 0)

    k = 0
    for i in eachindex(hx)
        if hx[i] != 0
            push!(ops, (-1, i))
            push!(p, h)
        end
    end

    # only take J_ij terms from upper diagonal since it's symmetric
    for j in axes(J, 2), i in axes(J, 1)
        if i < j
            if J[i, j] != 0
                push!(ops, op)
                push!(p, 2*abs(J))
            end
        end
    end

    return ops, p
end

function make_prob_vector(bond_spins::Vector{NTuple{2,Int}}, Ns::Int, J, h)
    ops = Vector{NTuple{2, Int}}(undef, 0)
    p = Vector{Float64}(undef, 0)

    if !iszero(h)
        for i in 1:Ns
            push!(ops, (-1, i))
            push!(p, h)
        end
    end

    if !iszero(J)
        for op in bond_spins
            push!(ops, op)
            push!(p, 2*abs(J))
        end
    end

    return ops, p
end


struct OperatorSampler{N, T <: Real}
    operators::Vector{NTuple{N, Int}}
    prob_heap::Vector{T}
    op_indices::Dict{NTuple{N, Int}, Int}
end

function OperatorSampler(operators::Vector{NTuple{N, Int}}, p::AbstractVector{T}) where {N, T <: Real}
    L = length(p)
    d = 2 ^ ceil(Int, log2(L))

    heap = zeros(T, 2*d)
    heap[d : d + L - 1] = p

    @simd for i in (d-1):-1:1
        heap[i] = heap[2*i] + heap[2*i + 1]
    end

    op_indices = Dict(op => i for (i, op) in enumerate(operators))

    return OperatorSampler{N, T}(operators, heap, op_indices)
end

length(os::OperatorSampler) = length(os.operators)

function setindex!(os::OperatorSampler, v::Float64, i::Int)
    heap = os.prob_heap
    l = length(heap)
    j = (l÷2 - 1) + i
    heap[j] = v
    while j > 1
        j ÷= 2
        heap[j] = heap[2*j] + heap[2*j + 1]
    end
end

function getindex(os::OperatorSampler, i::Int)
    denom = os.prob_heap[1]
    return os.prob_heap[(length(os.prob_heap) ÷ 2 - 1) + i] / denom
end
getindex(os::OperatorSampler{N}, op::NTuple{N, Int}) where N = os[os.op_indices[op]]

function getindex(os::OperatorSampler, r::AbstractArray{Int})
    denom = os.prob_heap[1]
    d = (length(os.prob_heap) ÷ 2) - 1
    return os.prob_heap[d .+ r] / denom
end
getindex(os::OperatorSampler{N}, ops::Vector{NTuple{N, Int}}) where N = os[os.op_indices[ops]]


firstindex(::OperatorSampler) = 1
lastindex(os::OperatorSampler) = length(os)

function rand(os::OperatorSampler{N})::NTuple{N, Int} where N
    heap = os.prob_heap
    l = length(heap) ÷ 2
    r = rand() * heap[1]

    i = 1
    while i < l
        i *= 2
        left = heap[i]
        if r > left
            r -= left
            i += 1
        end
    end
    return os.operators[(i - l + 1)]
end


# TODO: hierarchical OperatorSampler
# procedure is essentially:
# - combine all equal probabilities into one entry (add them together)
# - make a dictionary, mapping the reduced probability vector's indices to lists
#     of equiprobable operators
# - sample from the probability vector using heap sampling, then uniformly
#     sample from the selected list of operators
# - useful if there's a lot of equal probabilities

# won't be able to easily update the probabilities (but as it stands, I don't
# think we're need that functionality)
