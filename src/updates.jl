# updates.jl
#
# Defines the functions that perform the diagonal update, and also
# that build the linked list and operator cluster update

#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j)
@inline isdiagonal(op::NTuple{2,Int}) = (op[1] != -2)
@inline isidentity(op::NTuple{2,Int}) = (op[1] == 0)
@inline issiteoperator(op::NTuple{2,Int}) = (op[1] < 0)
@inline isbondoperator(op::NTuple{2,Int}) = (op[1] > 0)


function mc_step!(f::Function, qmc_state::BinaryQMCState, H::TFIM)
    diagonal_update!(qmc_state, H)

    cluster_data = linked_list_update(qmc_state, H)

    f(cluster_data, qmc_state, H)

    cluster_update!(cluster_data, qmc_state, H)
end

mc_step!(qmc_state, H) = mc_step!((args...) -> nothing, qmc_state, H)

########################## finite-beta #######################################

function mc_step_beta!(f::Function, qmc_state::BinaryQMCState, H::TFIM, beta::Real; eq::Bool = false)
    num_ops = diagonal_update_beta!(qmc_state, H, beta; eq = eq)

    cluster_data = linked_list_update_beta(qmc_state, H)

    f(cluster_data, qmc_state, H)

    cluster_update_beta!(cluster_data, qmc_state, H)

    return num_ops
end

mc_step_beta!(qmc_state, H, beta; eq = false) = mc_step_beta!((args...) -> nothing, qmc_state, H, beta; eq = eq)


# returns true if operator insertion succeeded
function insert_diagonal_operator!(qmc_state::BinaryQMCState, H::TFIM, spin_prop, n)
    # rr = rand()
    # if rr < H.P_h  # probability to choose a single-site operator
    #     qmc_state.operator_list[n] = (-1, rand(1:H.Ns))
    #     return true
    # else
    #     site1, site2 = H.bond_spin[rand(1:H.Nb)]
    #     # spins at each end of the bond must be the same
    #     if spin_prop[site1] == spin_prop[site2]
    #         qmc_state.operator_list[n] = (site1, site2)
    #         return true
    #     end
    # end
    site1, site2 = rand(H.op_sampler)
    if site1 < 0 || spin_prop[site1] == spin_prop[site2]
        qmc_state.operator_list[n] = (site1, site2)
        return true
    else
        return false
    end
end


#############################################################################

function diagonal_update!(qmc_state::BinaryQMCState, H::TFIM)
    spin_prop = copy(qmc_state.left_config)  # the propagated spin state

    for (n, op) in enumerate(qmc_state.operator_list)
        if !isdiagonal(op)
            spin_prop[op[2]] ⊻= 1  # spinflip
        else
            success = false
            while !success
                success = insert_diagonal_operator!(qmc_state, H, spin_prop, n)
            end
        end
    end

    # DEBUG
    # if spin_prop != qmc_state.right_config  # check the spin propagation for error
    #     error("Basis state propagation error in diagonal update!")
    # end
end


#############################################################################

nullt = (0, 0, 0)  # a null tuple

function linked_list_update(qmc_state::BinaryQMCState, H::TFIM)
    Ns = nspins(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config

    len = 2 * Ns
    for op in qmc_state.operator_list
        if issiteoperator(op)
            len += 2
        else
            len += 4
        end
    end

    # initialize linked list data structures
    LinkList = zeros(Int, len)  # needed for cluster update
    LegType = falses(len)

    # A diagonal bond operator has non trivial associates for cluster building
    Associates = [nullt for _ in 1:len]

    First = collect(1:Ns)
    idx = 0

    # The first N elements of the linked list are the spins of the LHS basis state
    for i in 1:Ns
        idx += 1
        LegType[idx] = spin_left[i]
    end

    spin_prop = copy(spin_left)  # the propagated spin state

    # Now, add the 2M operators to the linked list. Each has either 2 or 4 legs
    @inbounds for op in qmc_state.operator_list
        if issiteoperator(op)
            site = op[2]
            # lower or left leg
            idx += 1
            LinkList[idx] = First[site]
            LegType[idx] = spin_prop[site]
            current_link = idx

            if !isdiagonal(op)  # off-diagonal site operator
                spin_prop[site] ⊻= 1  # spinflip
            end

            LinkList[First[site]] = current_link  # completes backwards link
            First[site] = current_link + 1

            # upper or right leg
            idx += 1
            LegType[idx] = spin_prop[site]
        else  # diagonal bond operator
            site1, site2 = op

            # lower left
            idx += 1
            LinkList[idx] = First[site1]
            LegType[idx] = spin_prop[site1]
            current_link = idx

            LinkList[First[site1]] = current_link  # completes backwards link
            First[site1] = current_link + 2
            vertex1 = current_link
            Associates[idx] = (vertex1 + 1, vertex1 + 2, vertex1 + 3)

            # lower right
            idx += 1
            LinkList[idx] = First[site2]
            LegType[idx] = spin_prop[site2]
            current_link = idx

            LinkList[First[site2]] = current_link  # completes backwards link
            First[site2] = current_link + 2
            Associates[idx] = (vertex1, vertex1 + 2, vertex1 + 3)

            # upper left
            idx += 1
            LegType[idx] = spin_prop[site1]
            Associates[idx] = (vertex1, vertex1 + 1, vertex1 + 3)

            # upper right
            idx += 1
            LegType[idx] = spin_prop[site2]
            Associates[idx] = (vertex1, vertex1 + 1, vertex1 + 2)
        end
    end

    # The last N elements of the linked list are the final spin state
    for i in 1:Ns
        idx += 1
        LinkList[idx] = First[i]
        LegType[idx] = spin_prop[i]
        LinkList[First[i]] = idx
    end

    # DEBUG
    # if spin_prop != spin_right
    #     @debug "Basis state propagation error: LINKED LIST"
    # end

    return ClusterData(LinkList, LegType, Associates, First, nothing)

end

#############################################################################

function cluster_update!(cluster_data::ClusterData, qmc_state::BinaryQMCState, H::TFIM)
    Ns = nspins(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config
    operator_list = qmc_state.operator_list

    LinkList = cluster_data.linked_list
    LegType = cluster_data.leg_types
    Associates = cluster_data.associates

    lsize = length(LinkList)

    in_cluster = zeros(Int, lsize)
    cstack = Stack{Int}()  # This is the stack of vertices in a cluster
    ccount = 0  # cluster number counter

    @inbounds for i in 1:lsize
        # Add a new leg onto the cluster
        if (in_cluster[i] == 0 && Associates[i] === nullt)
            ccount += 1
            push!(cstack, i)
            in_cluster[i] = ccount

            flip = rand(Bool)  # flip a coin for the SW cluster flip
            if flip
                LegType[i] ⊻= 1  # spinflip
            end

            while !isempty(cstack)
                leg = LinkList[pop!(cstack)]

                if in_cluster[leg] == 0
                    in_cluster[leg] = ccount  # add the new leg and flip it
                    if flip
                        LegType[leg] ⊻= 1
                    end

                    # now check all associates and add to cluster
                    assoc = Associates[leg]  # a 3-tuple
                    if assoc !== nullt
                        for a in assoc
                            push!(cstack, a)
                            in_cluster[a] = ccount
                            if flip
                                LegType[a] ⊻= 1
                            end
                        end
                    end
                end

            end
        end
    end

    #println(in_cluster)

    # map back basis states and operator list
    for i in 1:Ns
        spin_left[i] = LegType[i]  # left basis state
        spin_right[i] = LegType[lsize-Ns+i]  # right basis state
    end

    ocount = Ns + 1  # next on is leg Ns + 1
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

end

#############################################################################
#############  FINITE BETA FUNCTIONS BELOW ##################################
#############################################################################

function diagonal_update_beta!(qmc_state::BinaryQMCState, H::TFIM, beta::Real; eq::Bool = false)

    # define the Metropolis probability as a constant
    # https://pitp.phas.ubc.ca/confs/sherbrooke2012/archives/Melko_SSEQMC.pdf
    # equation 1.42
    P_norm = beta * H.P_normalization

    num_ids = count(isidentity, qmc_state.operator_list)
    P_remove = (num_ids + 1) / P_norm
    P_accept = P_norm / num_ids

    spin_prop = copy(qmc_state.left_config)  # the propagated spin state

    for (n, op) in enumerate(qmc_state.operator_list)
        if !isdiagonal(op)
            spin_prop[op[2]] ⊻= 1  # spinflip
        elseif !isidentity(op)
            if rand() < P_remove
                qmc_state.operator_list[n] = (0, 0)
                num_ids += 1
                P_remove = (num_ids + 1) / P_norm
                P_accept = P_norm / num_ids
            end
        else
            if rand() < P_accept
                success = insert_diagonal_operator!(qmc_state, H, spin_prop, n)

                if success
                    # save one operation lol
                    P_remove = num_ids / P_norm
                    num_ids -= 1
                    P_accept = P_norm / num_ids
                end
            end
        end
    end

    # DEBUG
    # if spin_prop != qmc_state.right_config  # check the spin propagation for error
    #     error("Basis state propagation error in diagonal update!")
    # end

    total_list_size = length(qmc_state.operator_list)
    num_ops = total_list_size - num_ids

    if eq && 1.2*num_ops > length(qmc_state.operator_list)
        resize_op_list!(qmc_state.operator_list, round(Int, 1.5*num_ops))
    end
    return num_ops
end

#############################################################################

function linked_list_update_beta(qmc_state::BinaryQMCState, H::TFIM)
    Ns = nspins(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config

    len = 0
    for op in qmc_state.operator_list
        if issiteoperator(op)
            len += 2
        elseif isbondoperator(op)
            len += 4
        end
    end

    # initialize linked list data structures
    LinkList = zeros(Int, len)  # needed for cluster update
    LegType = falses(len)

    # A diagonal bond operator has non trivial associates for cluster building
    Associates = [nullt for _ in 1:len]

    First = zeros(Int,Ns)  #initialize the First list
    Last = zeros(Int,Ns)   #initialize the Last list
    idx = 0

    spin_prop = copy(spin_left)  # the propagated spin state

    # Now, add the 2M operators to the linked list. Each has either 2 or 4 legs
    @inbounds for op in qmc_state.operator_list
        if issiteoperator(op)
            site = op[2]
            # lower or left leg
            idx += 1
            LinkList[idx] = First[site]
            LegType[idx] = spin_prop[site]
            current_link = idx

            if !isdiagonal(op)  # off-diagonal site operator
                spin_prop[site] ⊻= 1  # spinflip
            end

            if First[site] != 0
                LinkList[First[site]] = current_link  # completes backwards link
            else
                Last[site] = current_link
            end
            First[site] = current_link + 1

            # upper or right leg
            idx += 1
            LegType[idx] = spin_prop[site]
        elseif isbondoperator(op)  # diagonal bond operator
            site1, site2 = op

            # lower left
            idx += 1
            LinkList[idx] = First[site1]
            LegType[idx] = spin_prop[site1]
            current_link = idx

			if First[site1] != 0
                LinkList[First[site1]] = current_link  # completes backwards link
            else
                Last[site1] = current_link
            end

            First[site1] = current_link + 2
            vertex1 = current_link
            Associates[idx] = (vertex1 + 1, vertex1 + 2, vertex1 + 3)

            # lower right
            idx += 1
            LinkList[idx] = First[site2]
            LegType[idx] = spin_prop[site2]
            current_link = idx

			if First[site2] != 0
                LinkList[First[site2]] = current_link  # completes backwards link
            else
                Last[site2] = current_link
            end

            First[site2] = current_link + 2
            Associates[idx] = (vertex1, vertex1 + 2, vertex1 + 3)

            # upper left
            idx += 1
            LegType[idx] = spin_prop[site1]
            Associates[idx] = (vertex1, vertex1 + 1, vertex1 + 3)

            # upper right
            idx += 1
            LegType[idx] = spin_prop[site2]
            Associates[idx] = (vertex1, vertex1 + 1, vertex1 + 2)
        end
    end

    #Periodic boundary conditions for finite-beta
    for i in 1:Ns
		if First[i] != 0  #This might be encountered at high temperatures
            LinkList[First[i]] = Last[i]
            LinkList[Last[i]] = First[i]
        end
    end

    # DEBUG
     if spin_prop != spin_right
         @debug "Basis state propagation error: LINKED LIST"
     end

    return ClusterData(LinkList, LegType, Associates,First,Last)

end

#############################################################################

function cluster_update_beta!(cluster_data::ClusterData, qmc_state::BinaryQMCState, H::TFIM)
    Ns = nspins(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config
    operator_list = qmc_state.operator_list

    LinkList = cluster_data.linked_list
    LegType = cluster_data.leg_types
    Associates = cluster_data.associates

    lsize = length(LinkList)

    in_cluster = zeros(Int, lsize)
    cstack = Stack{Int}()  # This is the stack of vertices in a cluster
    ccount = 0  # cluster number counter

    @inbounds for i in 1:lsize
        # Add a new leg onto the cluster
        if (in_cluster[i] == 0 && Associates[i] === nullt)
            ccount += 1
            push!(cstack, i)
            in_cluster[i] = ccount

            flip = rand(Bool)  # flip a coin for the SW cluster flip
            if flip
                LegType[i] ⊻= 1  # spinflip
            end

            while !isempty(cstack)
                leg = LinkList[pop!(cstack)]

                if in_cluster[leg] == 0
                    in_cluster[leg] = ccount  # add the new leg and flip it
                    if flip
                        LegType[leg] ⊻= 1
                    end

                    # now check all associates and add to cluster
                    assoc = Associates[leg]  # a 3-tuple
                    if assoc !== nullt
                        for a in assoc
                            push!(cstack, a)
                            in_cluster[a] = ccount
                            if flip
                                LegType[a] ⊻= 1
                            end
                        end
                    end
                end

            end
        end
    end

    #println(in_cluster)

    # map back basis states and operator list
    First = cluster_data.first
    Last = cluster_data.last
    for i in 1:Ns
        if First[i] != 0
            spin_left[i] = LegType[Last[i]]  # left basis state
            spin_right[i] = LegType[First[i]]  # right basis state
        else
			spin_left[i] = rand(Bool)
			spin_right[i] = spin_left[i]   #randomly flip spins not connected to operators
		end

    end

    ocount = 1  # first leg
    @inbounds for (n, op) in enumerate(operator_list)
        if isbondoperator(op)
            ocount += 4
        elseif !isidentity(op)
            if LegType[ocount] == LegType[ocount+1]  # diagonal
                operator_list[n] = (-1, op[2])
            else  # off-diagonal
                operator_list[n] = (-2, op[2])
            end
            ocount += 2
        end
    end

	#println(operator_list)

end
