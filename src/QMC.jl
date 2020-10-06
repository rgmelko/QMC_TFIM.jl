module QMC

using ProgressMeter

using Measurements
using Statistics
using FFTW

using DelimitedFiles
using JLD2
using Printf

using DataStructures
using SparseArrays

import Base: zero
import Base: length, size, eltype, setindex!, getindex, firstindex, lastindex, rand



export BinaryQMCState, Hamiltonian, TFIM, nspins, nbonds, ClusterData, mc_step!, mc_step_beta!,
        sample, simulation_cell, magnetization, num_single_site_diag, num_single_site_offdiag,
        num_single_site, autocorrelation, correlation_time, jackknife, mean_and_stderr,
        lattice_bond_spins


include("lattice.jl")
include("op_sampler.jl")
include("hamiltonian.jl")
include("qmc_state.jl")
include("measurements.jl")
include("updates.jl")
include("error.jl")


end