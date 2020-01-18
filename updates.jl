# updates.jl
#
# Defines the functions that perform the diagonal update, and also
# that build the linked list and operator cluster update

lsize = 0
nullt = (0, 0, 0) #a null tuple

############################ FUNCTIONS ######################################

#Diagonal update
function diagonal_update!(operator_list, h_x, J_, nSpin, nBond, spin_left, spin_right)

    #define the Metropolis probability as a constant
    #https://pitp.phas.ubc.ca/confs/sherbrooke2012/archives/Melko_SSEQMC.pdf
    #equation 1.43
    P_h = h_x * nSpin / (h_x * nSpin + 2.0 * J_ * nBond) #J=1.0 tested only

    spin_prop = copy(spin_left)  #the propagated spin state
    for i in 1:2*M  #size of the operator list

        if operator_list[i, 1] == -2
            spin_prop[operator_list[i, 2]] ⊻= 1 #spinflip
        else
            while true
                rr = rand()
                if rr < P_h #probability to choose a single-site operator
                    operator_list[i, 1] = -1
                    operator_list[i, 2] = rand(1:nSpin)  # site index
                    break
                else
                    bond = rand(1:nBond)
                    # spins at each end of the bond must be the same
                    if spin_prop[bond_spin[bond, 1]] == spin_prop[bond_spin[bond, 2]]
                        operator_list[i, :] = bond_spin[bond, :]
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
function LinkedList()

    #initialize linked list data structures
    global LinkList = zeros(Int, 0)  #needed for cluster update
    global LegType = zeros(Int, 0)

    #A diagonal bond operator has non trivial associates for cluster building
    #Associates = zeros((Int,Int,Int),0)
    global Associates = []

    First = zeros(Int, 0)  #scope is this function only

    #The first N elements of the linked list are the spins of the LHS basis state
    for i in 1:nSpin
        push!(First, i)
        push!(LinkList, -99) #we don't yet know what these link to
        push!(LegType, spin_left[i])
        push!(Associates, nullt) #no nontrivial associates for train wf spins
    end #i

    spin_prop = copy(spin_left)  #the propagated spin state

    #Now, add the 2M operators to the linked list.  Each has either 2 or 4 legs
    for i in 1:2*M  #size of the operator list

        if operator_list[i, 1] == -2  #off-diagonal site operator
            site = operator_list[i, 2]
            #lower or left leg
            push!(LinkList, First[site])
            push!(LegType, spin_prop[site]) #the spin of the vertex leg
            spin_prop[site] = xor(spin_prop[site], 1) #spinflip
            current_link = size(LinkList)
            LinkList[First[site]] = current_link[1] #completes backwards link
            First[site] = current_link[1] + 1
            push!(Associates, nullt)
            #upper or right leg
            push!(LinkList, -99) #we don't yet know what this links to
            push!(LegType, spin_prop[site]) #the spin of the vertex leg
            push!(Associates, nullt)

        elseif operator_list[i, 1] == -1  #diagonal site operator
            site = operator_list[i, 2]
            #lower or left leg
            push!(LinkList, First[site])
            push!(LegType, spin_prop[site]) #the spin of the vertex leg
            current_link = size(LinkList)
            LinkList[First[site]] = current_link[1] #completes backwards link
            First[site] = current_link[1] + 1
            push!(Associates, nullt)
            #upper or right leg
            push!(LinkList, -99) #we don't yet know what this links to
            push!(LegType, spin_prop[site]) #the spin of the vertex leg
            push!(Associates, nullt)

        else  #diagonal bond operator
            #lower left
            site1 = operator_list[i, 1]
            push!(LinkList, First[site1])
            push!(LegType, spin_prop[site1]) #the spin of the vertex leg
            current_link = size(LinkList)
            LinkList[First[site1]] = current_link[1] #completes backwards link
            First[site1] = current_link[1] + 2
            vertex1 = current_link[1]
            push!(Associates, (vertex1 + 1, vertex1 + 2, vertex1 + 3))
            #lower right
            site2 = operator_list[i, 2]
            push!(LinkList, First[site2])
            push!(LegType, spin_prop[site2]) #the spin of the vertex leg
            current_link = size(LinkList)
            LinkList[First[site2]] = current_link[1] #completes backwards link
            First[site2] = current_link[1] + 2
            push!(Associates, (vertex1, vertex1 + 2, vertex1 + 3))
            #upper left
            push!(LinkList, -99) #we don't yet know what this links to
            push!(LegType, spin_prop[site1]) #the spin of the vertex leg
            push!(Associates, (vertex1, vertex1 + 1, vertex1 + 3))
            #upper right
            push!(LinkList, -99) #we don't yet know what this links to
            push!(LegType, spin_prop[site2]) #the spin of the vertex leg
            push!(Associates, (vertex1, vertex1 + 1, vertex1 + 2))

        end #if

    end #i

    #The last N elements of the linked list are the final spin state
    for i in 1:nSpin
        push!(LinkList, First[i])
        push!(LegType, spin_prop[i])
        current_link = size(LinkList)
        LinkList[First[i]] = current_link[1]
        push!(Associates, nullt)
    end #i

    global lsize = size(LinkList)

    #DEBUG
    if spin_prop != spin_right
        println("Basis state propagation error: LINKED LIST")
    end

end #LinkedList

#############################################################################

#ClusterUpdate
function ClusterUpdate()

    #lsize is the size of the linked list
    in_cluster = zeros(Int, lsize[1])

    cstack = zeros(Int, 0)  #This is the stack of vertices in a cluster

    ccount = 0 #cluster number counter
    for i in 1:lsize[1]

        #Add a new leg onto the cluster
        if (in_cluster[i] == 0 && Associates[i] == nullt)

            ccount += 1
            push!(cstack, i)
            in_cluster[cstack[end]] = ccount

            flip = rand(Bool) #flip a coin for the SW cluster flip
            if flip == true
                LegType[cstack[end]] ⊻= 1 #spinflip
            end

            while isempty(cstack) == false

                leg = LinkList[cstack[end]]
                pop!(cstack)

                if in_cluster[leg] == 0

                    in_cluster[leg] = ccount #add the new leg and flip it
                    if flip == true
                        LegType[leg] ⊻= 1
                    end
                    #now check all associates and add to cluster
                    assoc = Associates[leg] #a 3-element array
                    if assoc != nullt
                        push!(cstack, assoc...)
                        in_cluster[collect(assoc)] .= ccount
                        if flip == true
                            LegType[collect(assoc)] .⊻= 1
                        end
                    end #if

                end #if in_cluster == 0
            end #while
        end #if
    end #for i

    #DEBUG
    #for i = 1:lsize[1]
    #    println(i," ",LegType[i]," ",in_cluster[i])
    #end

    #map back basis states and operator list
    ocount = 0
    for i in 1:nSpin
        spin_left[i] = LegType[i]  #left basis state
        ocount += 1
    end

    ocount += 1  #next on is leg nSpin + 1
    for i in 1:2*M
        if operator_list[i, 1] != -2 && operator_list[i, 1] != -1
            ocount += 4
        else
            if LegType[ocount] == LegType[ocount+1]  #diagonal
                operator_list[i, 1] = -1
            else
                operator_list[i, 1] = -2 #off-diagonal
            end
            ocount += 2
        end
    end

    for i in 1:nSpin
        spin_right[i] = LegType[lsize[1]-nSpin+i]  #left basis state
    end


end #ClusterUpdate
