# main.jl
#
# A projector QMC program for the TFIM

using Random
using ProgressMeter
Random.seed!(1234)

using DelimitedFiles
using Lattices

include("measurements.jl") # functions for the Monte Carlo updates


Dim = 1
nX = 6
PBC = true

BC = PBC ? Periodic : Fixed
BC_name = PBC ? "PBC" : "OBC"

if Dim == 1
    lattice = Chain(nX; boundary = BC)
else
    lattice = Square(nX, nX; boundary = BC)
end

# Projector parameters
M = 500 # length of the projector operator_list is 2M
h = 1.0
J_ = 1.0
MCS = 50000 # the number of Monte Carlo steps

root = "./data/$(Dim)D/$(nX)/$(BC_name)/h$(h)/"
mkpath(root)
output_file = "$(root)samples.txt"

# *******  Globals
include("updates.jl") # functions for the Monte Carlo updates

H = TFIM(lattice, h, J_)
qmc_state = BinaryQMCState(H, M)

@showprogress "warming up..." for i in 1:div(MCS, 5)  # Equilibration
    mc_step!(qmc_state, H)
end

measurements = zeros(Int, MCS, nspins(H))

mags = zeros(MCS)
ns = zeros(MCS)

@showprogress "MCMC..." for i in 1:MCS # Monte Carlo Production Steps
    mc_step!(qmc_state, H)

    spin_prop = sample(qmc_state.left_config, qmc_state.operator_list)
    measurements[i, :] = spin_prop

    ns[i] = num_single_site_diag(qmc_state.operator_list)
    mags[i] = magnetization(spin_prop)
end

open(output_file, "w") do io
    writedlm(io, measurements, " ")
end

n_inv = 1. / mean(ns)

@info "mean M is:\t" mean(mags)
@info "mean |M| is:\t" mean(abs, mags)
@info "mean M²:\t"  mean(x->x^2, mags)
@info "Energy is:\t" h*(n_inv - 1) - J_*nbonds(H)/nspins(H)
