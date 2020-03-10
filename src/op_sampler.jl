# (-3, i) = long field
# (-2, i) = transverse field off-diag
# (-1, i) = transverse field diag
# (0, 0) = id
# (i, j) = diag bond op


function make_prob_vector(J::AbstractMatrix, hx::AbstractVector)
    ops = repeat([(0, 0)], length(J) + length(hx))
    p = zeros(length(ops))

    k = 0
    for i in eachindex(hx)
        k += 1
        if hx[i] != 0
            ops[k] = (-1, i)
            p[k] = hx[i]
        end
    end

    # only take J_ij terms from upper diagonal since it's symmetric
    for j in axes(J, 2), i in axes(J, 1)
        if i < j
            k += 1
            if J[i, j] != 0
                ops[k] = (i, j)
                p[k] = 2*abs(J[i, j])
            end
        end
    end

    return ops, p
end

struct OperatorSampler{N}
    operators::Vector{NTuple{N, Int}}
    prob_heap::Vector{Float64}
end

function OperatorSampler(operators::Vector{NTuple{N, Int}}, p::AbstractVector{T}) where {N, T <: Real}
    L = length(p)
    d = 2 ^ ceil(Int, log2(L))

    heap = zeros(T, 2*d)
    heap[d : d + L - 1] = p

    @simd for i in (d-1):-1:1
        heap[i] = heap[2*i] + heap[2*i + 1]
    end

    return OperatorSampler{N}(operators, heap)
end

function rand(os::OperatorSampler)
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
    return os.operators[i - l + 1]
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
