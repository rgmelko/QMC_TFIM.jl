# measurements.jl
#
# Defines estimators and provides measurements
using Statistics


function sample(spin_left, operator_list)
    M = length(operator_list) ÷ 2
    spin_prop = copy(spin_left)

    for o in operator_list[1:M] #propagate half the list only (to the middle)
        if issiteoperator(o) && !isdiagonal(o)
            spin_prop[o.i] ⊻= 1 #spinflip
        end
    end
    return spin_prop
end

magnetization(spin_prop) = mean(x->2x - 1, spin_prop)

function interaction(J, spin_prop, BC::Type{Periodic})
    nSpin = length(spin_prop)
    spins = (2 * spin_prop) .- 1
    return -J*sum(spins .* spins[(2:nSpin+1) .% nSpin])
end

function interaction(J, spin_prop, BC::Type{Fixed})
    nSpin = length(spin_prop)
    spins = (2 * spin_prop) .- 1
    return -J*sum(spins[1:nSpin-1] .* spins[2:nSpin])
end

function interaction_ops(J, operator_list)
    return -J*mean(isbondoperator, operator_list)*nSpin
end


function num_single_site_diag(operator_list)
    return mean(x -> issiteoperator(x) && isdiagonal(x), operator_list)
end


function num_single_site_offdiag(operator_list)
    return mean(x -> issiteoperator(x) && !isdiagonal(x), operator_list)
end


function num_single_site(operator_list)
    return mean(issiteoperator, operator_list)
end
