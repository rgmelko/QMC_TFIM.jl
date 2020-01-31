# main.jl
#
# A projector QMC program for the TFIM

using Random
using ProgressMeter
Random.seed!(1234)

using DelimitedFiles
using JLD2

using Lattices

Dim = 1
nX = 6
PBC = true
h = 1.5
J_ = 2.0

BC = PBC ? Periodic : Fixed
BC_name = PBC ? "PBC" : "OBC"

if Dim == 1
    lattice = Chain(nX; boundary = BC)
else
    lattice = Square(nX, nX; boundary = BC)
end

# MC parameters
M = 700 # length of the projector operator_list is 2M
MCS = 100000 # the number of Monte Carlo steps
EQ_MCS = div(MCS, 10)

root = "./data/$(Dim)D/$(nX)/$(BC_name)/h$(h)/"
mkpath(root)
samples_file = "$(root)samples.txt"
qmc_state_file = "$(root)state.jld2"

# *******  Globals
include("updates.jl") # functions for the Monte Carlo updates

H = TFIM(lattice, h, J_)
qmc_state = BinaryQMCState(H, M)

@showprogress "warming up..." for i in 1:EQ_MCS  # Equilibration
    mc_step!(qmc_state, H)
end

measurements = zeros(Int, MCS, nspins(H))
mags = zeros(MCS)
ns = zeros(MCS)

@showprogress "MCMC..." for i in 1:MCS # Monte Carlo Production Steps
    mc_step!(qmc_state, H) do cluster_data, qmc_state, H
        spin_prop = sample(qmc_state)
        measurements[i, :] = spin_prop

        ns[i] = num_single_site_diag(qmc_state.operator_list)
        mags[i] = magnetization(spin_prop)
    end
end

open(samples_file, "w") do io
    writedlm(io, measurements, " ")
end

@save qmc_state_file qmc_state

n_inv = 1.0 / mean(ns)

@info "mean M is:" mean(mags)
@info "mean |M| is:" mean(abs, mags)
@info "mean M²:" mean(x -> x^2, mags)
@info "Energy is:" J_ * nbonds(H) / nspins(H) - h * (n_inv - 1)

@info "Operator list length:" 2*M
@info "Number of MC steps:" MCS
@info "Number of equilibration steps:" EQ_MCS

@show Dim, nX, PBC, h, J_

@info "Samples outputted to file:" samples_file