# main.jl
#
# A projector QMC program for the TFIM

using Random
using ProgressMeter
Random.seed!(1234)

using DelimitedFiles

include("lattice.jl") # define the spatial lattice
include("updates.jl") # functions for the Monte Carlo updates
include("measurements.jl") # functions for the Monte Carlo updates


Dim = 1
nX = 4
PBC = false

BC = PBC ? Periodic : Fixed
BC_name = PBC ? "PBC" : "OBC"

# lattice = HyperRect(lattice_shape; boundary = BC)
if Dim == 1
    lattice = Chain(nX; boundary = BC)
else
    lattice = Square(nX, nX; boundary = BC)
end
bond_spin = lattice_bond_spins(lattice)

nSpin = prod(size(lattice))
nBond = size(bond_spin, 1)

# Projector parameters
M = 400 # length of the projector operator_list is 2M
h = 2.0
J_ = 1.0
MCS = 8000 # the number of Monte Carlo steps

root = "./data/$(Dim)D/$(nX)/$(BC_name)/h$(h)/"
mkpath(root)
output_file = "$(root)samples.txt"

# *******  Globals
spin_left = falses(nSpin) # left and right trail spin state
spin_right = falses(nSpin)
LinkList = []
LegType = []
Associates = []

function make_op_list(M)
    operator_list::Vector{SSEOperator} = [IdOperator(1, lattice) for _ in 1:2M]
    return operator_list
end

operator_list = make_op_list(M)
#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j + 1)
# *******  Globals

@showprogress "warming up..." for i in 1:div(MCS, 5)  # Equilibration
    diagonal_update!(operator_list, lattice, h, J_, nSpin, nBond, spin_left, spin_right)
    LinkedList()
    ClusterUpdate(operator_list, lattice, nSpin, spin_left, spin_right, Associates, LinkList, LegType)
end

measurements = falses(MCS, nSpin)

mags = zeros(MCS)
energies = zeros(MCS)
E1 = zeros(MCS)
E2 = zeros(MCS)
inv = zeros(MCS)
ns = zeros(MCS)

@showprogress "MCMC..." for i in 1:MCS # Monte Carlo Production Steps
    # for j in 1:5
    diagonal_update!(operator_list, lattice, h, J_, nSpin, nBond, spin_left, spin_right)
    LinkedList()
    spin_prop = sample(spin_left, operator_list)
    ClusterUpdate(operator_list, lattice, nSpin, spin_left, spin_right, Associates, LinkList, LegType)

    measurements[i, :] = spin_prop

    E1[i], E2[i], inv[i], ns[i] = energy_abs_zero(h, J_, spin_prop, operator_list)
    mags[i] = magnetization(spin_prop)

    # end
    
end

Jenergy = mean(E1)
henergy = mean(E2)  # - J_ * nBond / nSpin
n1 = mean(inv)

@info "mean M is:\t" mean(mags)
@info "mean M²:\t"  mean(x->x^2, mags)
@info "Energy is:\t" Jenergy
@info "Energy is:\t" henergy
@info "1/n is:\t" n1

@info "Energy is:\t" henergy + Jenergy - h - J_ * nBond / nSpin


@info J_ * nBond / nSpin
@info h