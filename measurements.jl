# measurements.jl
#
# Defines estimators and provides measurements
using Statistics
using FFTW

include("qmc.jl")

function sample(qmc_state::BinaryQMCState)
    operator_list = qmc_state.operator_list

    M = length(operator_list) ÷ 2
    spin_prop = copy(qmc_state.left_config)

    for op in operator_list[1:M] #propagate half the list only (to the middle)
        if issiteoperator(op) && !isdiagonal(op)
            spin_prop[op[2]] ⊻= 1 #spinflip
        end
    end
    return spin_prop
end

magnetization(spin_prop) = mean(x -> 2x - 1, spin_prop)

num_single_site_diag(operator_list) = mean(x -> issiteoperator(x) && isdiagonal(x), operator_list)
num_single_site_offdiag(operator_list) = mean(x -> issiteoperator(x) && !isdiagonal(x), operator_list)
num_single_site(operator_list) = mean(issiteoperator, operator_list)


function autocorrelation(m::Vector)
    N = length(m)

    m′ = m .- mean(m)
    m′ = vcat(m′, zeros(N))
    mw = fft(m′)
    s = abs2.(mw)

    chi = real(ifft(s)[1:N])

    for i in 1:N
        chi[i] /= (2*N)  # normalize FFT
        chi[i] /= (N - i - 1)
    end
    return chi
end


# use method explained by Sokal to estimate correlation time
# https://pdfs.semanticscholar.org/0bfe/9e3db30605fe2d4d26e1a288a5e2997e7225.pdf
function correlation_time(m::Vector)
    ac = autocorrelation(m)
    ac_0 = ac[1]

    corr_time = 1.0
    @inbounds for M in axes(ac, 1)
        corr_time += 2*(ac[M] / ac_0)
        if M >= 5*corr_time
            break
        end
    end

    return corr_time
end