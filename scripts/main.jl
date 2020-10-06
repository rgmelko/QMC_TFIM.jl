using DrWatson
@quickactivate "QMC"

using QMC

# main.jl
#
# A projector QMC program for the TFIM

using Random
Random.seed!(1234)

using ProgressMeter

using Measurements
using Statistics
using FFTW

using DelimitedFiles
using JLD2
using Printf
using Lattices

using DataStructures
using ArgParse


###############################################################################

function init_mc_cli(parsed_args)
    PBC = parsed_args["periodic"]
    h = parsed_args["field"]
    J = parsed_args["interaction"]

    Dim = length(parsed_args["dims"])
    nX = parsed_args["dims"]

    BC = PBC ? Periodic : Fixed
    BC_name = PBC ? "PBC" : "OBC"

    if Dim == 1
        nX = nX[1]
        lattice = Chain(nX; boundary = BC)
    elseif Dim == 2
        lattice = Square(nX...; boundary = BC)
    else
        error("Unsupported number of dimensions")
    end

    # MC parameters
    M = parsed_args["M"] # length of the operator_list is 2M
    MCS = parsed_args["measurements"] # the number of samples to record
    EQ_MCS = div(MCS, 10)
    skip = parsed_args["skip"]  # number of MC steps to perform between each msmt

    # path = "$(Dim)D/$(nX)/$(BC_name)/J$(J)/h$(h)/skip$(skip)/"
    d = @ntuple Dim nX BC_name J h skip M
    mc_opts = @ntuple M MCS EQ_MCS skip
    bond_spin = lattice_bond_spins(lattice)

    Ns, Nb = length(lattice), length(bond_spin)
    H = TFIM(bond_spin, Ns, Nb, h, J)
    qmc_state = BinaryQMCState(H, M)

    return H, qmc_state, savename(d; digits = 4), mc_opts
end


function make_info_file(info_file, samples_file, mc_opts, observables, corr_time)
    M, MCS, EQ_MCS, skip = mc_opts
    mag, abs_mag, mag_sqr, energy = observables

    open(info_file, "w") do file_io
        streams = [Base.stdout, file_io]

        for io in streams
            @printf(io, "⟨M⟩   = % .16f +/- %.16f\n", mag.val, mag.err)
            @printf(io, "⟨|M|⟩ = % .16f +/- %.16f\n", abs_mag.val, abs_mag.err)
            @printf(io, "⟨M^2⟩ = % .16f +/- %.16f\n", mag_sqr.val, mag_sqr.err)
            @printf(io, "⟨E⟩   = % .16f +/- %.16f\n\n", energy.val, energy.err)

            println(io, "Correlation time: $(corr_time)\n")

            println(io, "Operator list length: $(2*M)")
            println(io, "Number of MC measurements: $(MCS)")
            println(io, "Number of equilibration steps: $(EQ_MCS)")
            println(io, "Number of skips between measurements: $(skip)\n")

            println(io, "Samples outputted to file: $(samples_file)")
        end
    end
end


function save_data(path, mc_opts, qmc_state, measurements, observables, corr_time)
    info_file = path * "_info.txt"
    samples_file = path * "_samples.txt"
    qmc_state_file = path * "_state.jld2"

    open(samples_file, "w") do io
        writedlm(samples_file, measurements, " ")
    end

    @time @save qmc_state_file qmc_state

    make_info_file(info_file, samples_file, mc_opts, observables, corr_time)
end


function mixedstate(parsed_args)
    H, qmc_state, sname, mc_opts = init_mc_cli(parsed_args)
    beta = parsed_args["beta"]

    M, MCS, EQ_MCS, skip = mc_opts

    mkpath(datadir("sims", "mixedstate"))
    path = datadir("sims", "mixedstate", savename((@ntuple beta), sname; digits = 4))

    measurements = zeros(Int, MCS, nspins(H))
    mags = zeros(MCS)
    ns = zeros(MCS)

    @showprogress "Warm up..." for i in 1:EQ_MCS
        mc_step_beta!(qmc_state, H, beta; eq = true)
    end

    @showprogress "MCMC...   " for i in 1:MCS # Monte Carlo Steps
        ns[i] = mc_step_beta!(qmc_state, H, beta) do cluster_data, qmc_state, H
            spin_prop = qmc_state.left_config
            measurements[i, :] = spin_prop
            mags[i] = magnetization(spin_prop)
        end

        for _ in 1:skip
            mc_step_beta!(qmc_state, H, beta)
        end
    end

    mag = mean_and_stderr(mags)
    abs_mag = mean_and_stderr(abs, mags)
    mag_sqr = mean_and_stderr(abs2, mags)

    energy = mean_and_stderr(x -> -x/beta, ns) + H.J*nbonds(H) + H.h*nspins(H)
    energy /= nspins(H)

    observables = (mag, abs_mag, mag_sqr, energy)

    # measure correlation time from equilibriation samples
    @time corr_time = correlation_time(mags .^ 2)

    save_data(path, mc_opts, qmc_state, measurements, observables, corr_time)
end



function groundstate(parsed_args)
    H, qmc_state, sname, mc_opts = init_mc_cli(parsed_args)

    M, MCS, EQ_MCS, skip = mc_opts

    mkpath(datadir("sims", "groundstate"))
    path = datadir("sims", "groundstate", sname)

    measurements = zeros(Int, MCS, nspins(H))
    mags = zeros(MCS)
    ns = zeros(MCS)

    @showprogress "Warm up..." for i in 1:EQ_MCS
        mc_step!(qmc_state, H)
    end

    @showprogress "MCMC...   " for i in 1:MCS # Monte Carlo Production Steps
        mc_step!(qmc_state, H) do cluster_data, qmc_state, H
            spin_prop = sample(qmc_state)
            measurements[i, :] = spin_prop

            ns[i] = num_single_site_diag(qmc_state.operator_list)
            mags[i] = magnetization(spin_prop)
        end

        for _ in 1:skip
            mc_step!(qmc_state, H)
        end
    end

    mag = mean_and_stderr(mags)
    abs_mag = mean_and_stderr(abs, mags)
    mag_sqr = mean_and_stderr(abs2, mags)

    @time energy = jackknife(ns) do n
        if H.h != 0
            (-H.h * ((1.0 / n) - 1)) + H.J * (nbonds(H) / nspins(H))
        else
            H.J * (nbonds(H) / nspins(H))
        end
    end

    observables = (mag, abs_mag, mag_sqr, energy)

    # measure correlation time from equilibriation samples
    @time corr_time = correlation_time(mags .^ 2)

    save_data(path, mc_opts, qmc_state, measurements, observables, corr_time)
end


###############################################################################


s = ArgParseSettings()


@add_arg_table! s begin
    "groundstate"
        help = "Use Projector SSE to simulate the ground state"
        action = :command
    "mixedstate"
        help = "Use vanilla SSE to simulate the system at non-zero temperature"
        action = :command
end


@add_arg_table! s["groundstate"] begin
    "dims"
        help = "The dimensions of the square lattice"
        required = true
        arg_type = Int
        nargs = '+'
    "--periodic", "-p"
        help = "Periodic BCs"
        action = :store_true
    "--field"
        help = "Strength of the transverse field"
        arg_type = Float64
        default = 1.0
    "--interaction", "-J"
        help = "Strength of the interaction"
        arg_type = Float64
        default = 1.0

    "-M"
        help = "Half-size of the operator list"
        arg_type = Int
        default = 1000

    "--measurements", "-n"
        help = "Number of samples to record"
        arg_type = Int
        default = 100_000
    "--skip", "-s"
        help = "Number of MC steps to perform between each measurement"
        arg_type = Int
        default = 0
end

import_settings!(s["mixedstate"], s["groundstate"])

@add_arg_table! s["mixedstate"] begin
    "--beta"
        help = "The inverse-temperature parameter for the simulation"
        arg_type = Float64
        default = 10.0
end


parsed_args = parse_args(ARGS, s)

if parsed_args["%COMMAND%"] == "groundstate"
    groundstate(parsed_args["groundstate"])
else
    mixedstate(parsed_args["mixedstate"])
end
