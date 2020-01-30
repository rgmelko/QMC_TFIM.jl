# updates.jl
#
# Defines the functions that perform the diagonal update, and also
# that build the linked list and operator cluster update

include("qmc.jl")


struct ClusterData
    linked_list::Vector{Int}
    leg_types::Vector{Int}
    associates::Vector{NTuple{3, Int}}
end


function mc_step!(f::Function, qmc_state::BinaryQMCState, H::TFIM)
    diagonal_update!(qmc_state, H)
    cluster_data = linked_list_update(qmc_state, H)

    f(cluster_data, qmc_state, H)

    cluster_update!(cluster_data, qmc_state, H)
end

mc_step!(qmc_state, H) = mc_step!((args...) -> nothing, qmc_state, H)



#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j)
@inline isdiagonal(op::NTuple{2, Int}) = (op[1] != -2)
@inline issiteoperator(op::NTuple{2, Int}) = (op[1] <= 0)
@inline isbondoperator(op::NTuple{2, Int}) = (op[1] > 0)


#Diagonal update
function diagonal_update!(qmc_state::BinaryQMCState, H::TFIM)
    lattice, bond_spin = H.lattice, H.bond_spin
    Ns, Nb = nspins(H), nbonds(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config
    operator_list = qmc_state.operator_list

    #define the Metropolis probability as a constant
    #https://pitp.phas.ubc.ca/confs/sherbrooke2012/archives/Melko_SSEQMC.pdf
    #equation 1.43
    P_h = H.h * Ns / (H.h * Ns + 2.0 * H.J * Nb)

    spin_prop = copy(spin_left)  #the propagated spin state
    @inbounds for (n, op) in enumerate(operator_list)

        if issiteoperator(op) && !isdiagonal(op)
            spin_prop[op[2]] ⊻= 1 #spinflip
        else
            while true
                rr = rand()
                if rr < P_h #probability to choose a single-site operator
                    operator_list[n] = (-1, rand(1:Ns))
                    break
                else
                    bond = rand(1:Nb)
                    # spins at each end of the bond must be the same
                    site1, site2 = bond_spin[bond, :]
                    if spin_prop[site1] == spin_prop[site2]
                        operator_list[n] = (site1, site2)
                        break
                    end
                end #if rr < P_h
            end #while
        end # if
    end #for i

    #DEBUG
    if spin_prop != spin_right  #check the spin propagation for error
        error("Basis state propagation error in diagonal update!")
    end
end

#############################################################################

nullt = (0, 0, 0) #a null tuple

#LinkedList
function linked_list_update(qmc_state::BinaryQMCState, H::TFIM)
    Ns = nspins(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config
    operator_list = qmc_state.operator_list

    len = 2*Ns
    for op in operator_list
        if issiteoperator(op)
            len += 2
        else
            len += 4
        end
    end

    #initialize linked list data structures
    LinkList = zeros(Int, len)  #needed for cluster update
    LegType = zeros(Int, len)

    #A diagonal bond operator has non trivial associates for cluster building
    #Associates = zeros((Int,Int,Int),0)
    Associates = collect(repeat([nullt], len))

    First = collect(1:Ns)
    idx = 0

    #The first N elements of the linked list are the spins of the LHS basis state
    for i in 1:Ns
        idx += 1
        LegType[idx] = spin_left[i]
    end

    spin_prop = copy(spin_left)  #the propagated spin state

    #Now, add the 2M operators to the linked list.  Each has either 2 or 4 legs
    @inbounds for op in operator_list

        if issiteoperator(op)
            site = op[2]
            #lower or left leg
            idx += 1
            LinkList[idx] = First[site]
            LegType[idx] = spin_prop[site]
            current_link = idx

            if !isdiagonal(op)  #off-diagonal site operator
                spin_prop[site] ⊻= 1 #spinflip
            end

            LinkList[First[site]] = current_link #completes backwards link
            First[site] = current_link + 1

            #upper or right leg
            idx += 1
            LegType[idx] = spin_prop[site]
        else  #diagonal bond operator
            site1, site2 = op
            #lower left
            idx += 1
            LinkList[idx] = First[site1]
            LegType[idx] = spin_prop[site1]
            current_link = idx

            LinkList[First[site1]] = current_link #completes backwards link
            First[site1] = current_link + 2
            vertex1 = current_link
            Associates[idx] = (vertex1 + 1, vertex1 + 2, vertex1 + 3)

            #lower right
            idx += 1
            LinkList[idx] = First[site2]
            LegType[idx] = spin_prop[site2]
            current_link = idx

            LinkList[First[site2]] = current_link #completes backwards link
            First[site2] = current_link + 2
            Associates[idx] = (vertex1, vertex1 + 2, vertex1 + 3)

            #upper left
            idx += 1
            LegType[idx] = spin_prop[site1]
            Associates[idx] = (vertex1, vertex1 + 1, vertex1 + 3)

            #upper right
            idx += 1
            LegType[idx] = spin_prop[site2]
            Associates[idx] = (vertex1, vertex1 + 1, vertex1 + 2)
        end #if
    end #i

    #The last N elements of the linked list are the final spin state
    for i in 1:Ns
        idx += 1
        LinkList[idx] = First[i]
        LegType[idx] = spin_prop[i]
        LinkList[First[i]] = idx
    end

    #DEBUG
    if spin_prop != spin_right
        @debug "Basis state propagation error: LINKED LIST"
    end

    return ClusterData(LinkList, LegType, Associates)

end #LinkedList

#############################################################################

#ClusterUpdate
function cluster_update!(cluster_data::ClusterData, qmc_state::BinaryQMCState, H::TFIM)
    lattice = H.lattice
    Ns = nspins(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config
    operator_list = qmc_state.operator_list

    LinkList = cluster_data.linked_list
    LegType = cluster_data.leg_types
    Associates = cluster_data.associates

    lsize = length(LinkList)

    in_cluster = zeros(Int, lsize)
    cstack = zeros(Int, 0)  #This is the stack of vertices in a cluster
    ccount = 0 #cluster number counter
    @inbounds for i in 1:lsize
        #Add a new leg onto the cluster
        if (in_cluster[i] == 0 && Associates[i] == nullt)
            ccount += 1
            push!(cstack, i)
            in_cluster[i] = ccount

            flip = rand(Bool) #flip a coin for the SW cluster flip
            if flip
                LegType[i] ⊻= 1 #spinflip
            end
            while !isempty(cstack)

                leg = LinkList[pop!(cstack)]

                if in_cluster[leg] == 0
                    in_cluster[leg] = ccount #add the new leg and flip it
                    if flip
                        LegType[leg] ⊻= 1
                    end
                    #now check all associates and add to cluster
                    assoc = Associates[leg] #a 3-element array
                    if assoc != nullt
                        for a in assoc
                            push!(cstack, a)
                            in_cluster[a] = ccount
                            if flip
                                LegType[a] ⊻= 1
                            end
                        end
                    end #if
                end #if in_cluster == 0
            end #while
        end #if
    end #for i

    #map back basis states and operator list
    for i in 1:Ns
        spin_left[i] = LegType[i]  #left basis state
        spin_right[i] = LegType[lsize-Ns+i]  #right basis state
    end

    ocount = Ns + 1  #next on is leg nSpin + 1
    @inbounds for (n, op) in enumerate(operator_list)
        if isbondoperator(op)
            ocount += 4
        else
            if LegType[ocount] == LegType[ocount+1]  # diagonal
                operator_list[n] = (-1, op[2])
            else  # off-diagonal
                operator_list[n] = (-2, op[2])
            end
            ocount += 2
        end
    end

end #ClusterUpdate
