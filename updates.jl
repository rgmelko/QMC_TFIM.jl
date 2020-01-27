# updates.jl
#
# Defines the functions that perform the diagonal update, and also
# that build the linked list and operator cluster update

nullt = (0, 0, 0) #a null tuple

############################ FUNCTIONS ######################################
include("operators.jl")
include("qmc.jl")


struct ClusterData
    linked_list::Vector{Int}
    leg_types::Vector{Int}
    associates::Vector{NTuple{3, Int}}
end


@inline function mc_step!(qmc_state::BinaryQMCState, H::TFIM)
    diagonal_update!(qmc_state, H)
    cluster_data = linked_list_update(qmc_state, H)
    cluster_update!(cluster_data, qmc_state, H)
end

#Diagonal update
function diagonal_update!(qmc_state::BinaryQMCState, H::TFIM)
    lattice, bond_spin = H.lattice, H.bond_spin
    Ns, Nb = nspins(H), nbonds(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config
    operator_list = qmc_state.operator_list

    #define the Metropolis probability as a constant
    #https://pitp.phas.ubc.ca/confs/sherbrooke2012/archives/Melko_SSEQMC.pdf
    #equation 1.43
    P_h = H.h * Ns / (H.h * Ns + 2.0 * H.J * Nb) #J=1.0 tested only

    spin_prop = copy(spin_left)  #the propagated spin state
    for (n, op) in enumerate(operator_list)  #size of the operator list

        if issiteoperator(op) && !isdiagonal(op)
            spin_prop[op.i] ⊻= 1 #spinflip
        else
            while true
                rr = rand()
                if rr < P_h #probability to choose a single-site operator
                    operator_list[n] = SiteOperator{Diagonal}(rand(1:Ns), lattice)
                    break
                else
                    bond = rand(1:Nb)
                    # spins at each end of the bond must be the same
                    if spin_prop[bond_spin[bond, 1]] == spin_prop[bond_spin[bond, 2]]
                        operator_list[n] = BondOperator{Diagonal}(bond_spin[bond, :]..., lattice)
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

#LinkedList
function linked_list_update(qmc_state::BinaryQMCState, H::TFIM)
    Ns = nspins(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config
    operator_list = qmc_state.operator_list

    len = 2*Ns
    for op in operator_list
        if issiteoperator(op)
            len += 2
        elseif isbondoperator(op)
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
    for (n, op) in enumerate(operator_list)  #size of the operator list

        if issiteoperator(op)
            site = op.i
            #lower or left leg
            idx += 1
            LinkList[idx] = First[op.i]
            LegType[idx] = spin_prop[op.i]
            current_link = idx

            if !isdiagonal(op)  #off-diagonal site operator
                spin_prop[op.i] ⊻= 1 #spinflip
            end

            LinkList[First[op.i]] = current_link #completes backwards link
            First[op.i] = current_link + 1

            #upper or right leg
            idx += 1
            LegType[idx] = spin_prop[op.i]
        else  #diagonal bond operator
            #lower left
            idx += 1
            LinkList[idx] = First[op.i]
            LegType[idx] = spin_prop[op.i]
            current_link = idx

            LinkList[First[op.i]] = current_link #completes backwards link
            First[op.i] = current_link + 2
            vertex1 = current_link
            Associates[idx] = (vertex1 + 1, vertex1 + 2, vertex1 + 3)

            #lower right
            idx += 1
            LinkList[idx] = First[op.j]
            LegType[idx] = spin_prop[op.j]
            current_link = idx

            LinkList[First[op.j]] = current_link #completes backwards link
            First[op.j] = current_link + 2
            Associates[idx] = (vertex1, vertex1 + 2, vertex1 + 3)

            #upper left
            idx += 1
            LegType[idx] = spin_prop[op.i]
            Associates[idx] = (vertex1, vertex1 + 1, vertex1 + 3)

            #upper right
            idx += 1
            LegType[idx] = spin_prop[op.j]
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
    for i in 1:lsize
        #Add a new leg onto the cluster
        if (in_cluster[i] == 0 && Associates[i] == nullt)
            ccount += 1
            push!(cstack, i)
            in_cluster[cstack[end]] = ccount

            flip = rand(Bool) #flip a coin for the SW cluster flip
            if flip
                LegType[cstack[end]] ⊻= 1 #spinflip
            end
            while !isempty(cstack)

                leg = LinkList[cstack[end]]
                pop!(cstack)

                if in_cluster[leg] == 0
                    in_cluster[leg] = ccount #add the new leg and flip it
                    if flip
                        LegType[leg] ⊻= 1
                    end
                    #now check all associates and add to cluster
                    assoc = Associates[leg] #a 3-element array
                    if assoc != nullt
                        push!(cstack, assoc...)
                        assoc_l = collect(assoc)
                        in_cluster[assoc_l] .= ccount
                        if flip
                            LegType[assoc_l] .⊻= 1
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
    for (n, op) in enumerate(operator_list)
        if isbondoperator(op)
            ocount += 4
        else
            if LegType[ocount] == LegType[ocount+1]  #diagonal
                operator_list[n] = SiteOperator{Diagonal}(op.i, lattice)
            else
                operator_list[n] = SiteOperator{OffDiagonal}(op.i, lattice) #off-diagonal
            end
            ocount += 2
        end
    end

end #ClusterUpdate
