# main.jl
#
# A projector QMC program for the TFIM

using Random
using ProgressMeter
Random.seed!(1234)

using DelimitedFiles
using JLD2
using Printf

using Lattices

Dim = 1
nX = 6
PBC = false
h = 1.0
J_ = 1.0

BC = PBC ? Periodic : Fixed
BC_name = PBC ? "PBC" : "OBC"

if Dim == 1
    lattice = Chain(nX; boundary = BC)
else
    lattice = Square(nX, nX; boundary = BC)
end

# MC parameters
M = 200 # length of the projector operator_list is 2M
MCS = 300000 # the number of samples to record
EQ_MCS = div(MCS, 10)
skip = 0  # number of MC steps to perform between each msmt

root = "./data/$(Dim)D/$(nX)/$(BC_name)/J$(J_)/h$(h)/skip$(skip)/"
mkpath(root)
samples_file = "$(root)samples.txt"
info_file = "$(root)info.txt"
qmc_state_file = "$(root)state.jld2"

# *******  Globals
include("updates.jl") # functions for the Monte Carlo updates

H = TFIM(lattice, h, J_)
qmc_state = BinaryQMCState(H, M)

## FINITE-BETA
beta = 20.0
for i in 1:EQ_MCS  # Equilibration
    mc_step_beta!(qmc_state, H, beta) 
end

measurements = zeros(Int, MCS, nspins(H))
mags = zeros(MCS)
ns = zeros(MCS)

@showprogress "MCMC..." for i in 1:MCS # Monte Carlo Steps
    mc_step_beta!(qmc_state, H, beta) do cluster_data, qmc_state, H 
        #mags[i] = magnetization(qmc_state.left_config)	
        #spin_prop = sample(qmc_state)
        spin_prop = qmc_state.right_config
        mags[i] = magnetization(spin_prop)
    end
end

#include("error.jl")
#mag_sqr = mean_and_stderr(x -> x^2, mags)
#println(mag_sqr)

## ZERO-T PROJECTOR
#for i in 1:EQ_MCS  # Equilibration
#    mc_step!(qmc_state, H) 
#end
#
#measurements = zeros(Int, MCS, nspins(H))
#mags = zeros(MCS)
#ns = zeros(MCS)
#
#@showprogress "MCMC..." for i in 1:MCS # Monte Carlo Production Steps
#    mc_step!(qmc_state, H) do cluster_data, qmc_state, H
#        spin_prop = sample(qmc_state)
#        measurements[i, :] = spin_prop
#
#        ns[i] = num_single_site_diag(qmc_state.operator_list)
#        mags[i] = magnetization(spin_prop)
#    end
#
#    for _ in 1:skip
#        mc_step!(qmc_state, H)
#    end
#end

# open(samples_file, "w") do io
#     writedlm(io, measurements, " ")
# end
# 
# @time @save qmc_state_file qmc_state
# 
# energy(H::TFIM) = n -> (H.J * (nbonds(H) / nspins(H)) - H.h * ((1.0 / n) - 1))
 include("error.jl")
# 
# 
# mag = mean_and_stderr(mags)
# abs_mag = mean_and_stderr(abs, mags)
 mag_sqr = mean_and_stderr(x -> x^2, mags)
 println(mag_sqr)
 
# @time energy_ = jackknife(energy(H), ns)
# @time corr_time = correlation_time(mags)
# 
# println()
# 
# open(info_file, "w") do file_io
#     streams = [Base.stdout, file_io]
# 
#     for io in streams
#         @printf(io, "⟨M⟩   = % .16f +/- %.16f\n", mag.val, mag.err)
#         @printf(io, "⟨|M|⟩ = % .16f +/- %.16f\n", abs_mag.val, abs_mag.err)
#         @printf(io, "⟨M^2⟩ = % .16f +/- %.16f\n", mag_sqr.val, mag_sqr.err)
#         @printf(io, "⟨E⟩   = % .16f +/- %.16f\n\n", energy_.val, energy_.err)
# 
#         println(io, "Correlation time: $(corr_time)\n")
# 
#         println(io, "Operator list length: $(2*M)")
#         println(io, "Number of MC measurements: $(MCS)")
#         println(io, "Number of equilibration steps: $(EQ_MCS)")
#         println(io, "Number of skips between measurements: $(skip)\n")
# 
#         println(io, "Dim = $(Dim), nX = $(nX), PBC = $(PBC), h = $(h), J = $(J_)\n")
# 
#         println(io, "Samples outputted to file: $(samples_file)")
#     end
# end
