using Measurements

function mean_and_stderr(f::Function, x::Vector)
    y = f.(x)

    μ = mean(y)

    σ = std(y; mean=μ) / sqrt(length(x))

    return μ ± σ
end
mean_and_stderr(x::Vector) = mean_and_stderr(identity, x)


function jackknife(f::Function, x::Vector...)
    sum_x = [sum(x[i]) for i in 1:length(x)]
    N = length(x[1])

    f_J = zeros(N)
    @simd for i in 1:N
        x_J = [(sum_x[j] - x[j][i]) / (N - 1) for j in 1:length(x)]
        f_J[i] = f(x_J...)
    end

    μ = mean(f_J)
    σ = sqrt(N - 1) * std(f_J; mean=μ)

    return μ ± σ
end
jackknife(x::Vector) = jackknife(identity, x)