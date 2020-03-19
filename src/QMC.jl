module QMC

using ProgressMeter

using Measurements
using Statistics
using FFTW

using DelimitedFiles
using JLD2
using Printf
using Lattices

using DataStructures

import Base: zero


export BinaryQMCState, Hamiltonian, TFIM, nspins, nbonds, ClusterData, mc_step!, mc_step_beta!,
        sample, simultation_cell, magnetization, num_single_site_diag, num_single_site_offdiag,
        num_single_site, autocorrelation, correlation_time, jackknife, mean_and_stderr


include("lattice.jl")
include("hamiltonian.jl")
include("qmc_state.jl")
include("measurements.jl")
include("updates.jl")
include("error.jl")


end