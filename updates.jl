# updates.jl
#
# Defines the functions that perform the diagonal update, and also
# that build the linked list and operator cluster update

lsize = 0
nullt = (0, 0, 0) #a null tuple

############################ FUNCTIONS ######################################


abstract type OperatorForm end
struct Diagonal <: OperatorForm end
struct OffDiagonal <: OperatorForm end

abstract type SSEOperator{N, F <: OperatorForm, L <: BoundedLattice} end

struct SiteOperator{F, L} <: SSEOperator{1, F, L}
    i::Int
end
struct IdOperator{L} <: SSEOperator{1, Diagonal, L}
    i::Int
end
struct BondOperator{F, L} <: SSEOperator{2, F, L}
    i::Int
    j::Int
end

function SiteOperator{F}(i::Int, lat::L) where {L <: BoundedLattice, F <: OperatorForm}
    @assert 0 < i <= length(lat)
    return SiteOperator{F, L}(i)
end

function IdOperator(i::Int, lat::L) where L <: BoundedLattice
    @assert 0 < i <= length(lat)
    return IdOperator{L}(i)
end
function BondOperator{F}(i::Int, j::Int, lat::L) where {L <: BoundedLattice, F <: OperatorForm}
    @assert 0 < i <= length(lat)
    @assert 0 < j <= length(lat)
    return BondOperator{F, L}(i, j)
end

operatorform(::Type{<:SSEOperator{N, F}}) where {N, F} = F
operatorform(::T) where {T <: SSEOperator} = operatorform(T)

isdiagonal(T::Type{<:SSEOperator}) = (operatorform(T) <: Diagonal)
isdiagonal(::T) where {T <: SSEOperator} = isdiagonal(T)

issiteoperator(T::Type{<:SSEOperator}) = (T <: SiteOperator)
issiteoperator(::T) where {T <: SSEOperator} = issiteoperator(T)

isbondoperator(T::Type{<:SSEOperator}) = (T <: BondOperator)
isbondoperator(::T) where {T <: SSEOperator} = isbondoperator(T)

#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j + 1)




#Diagonal update
function diagonal_update!(operator_list, lattice, h, J_, nSpin, nBond, spin_left, spin_right)

    #define the Metropolis probability as a constant
    #https://pitp.phas.ubc.ca/confs/sherbrooke2012/archives/Melko_SSEQMC.pdf
    #equation 1.43
    P_h = h * nSpin / (h * nSpin + 2.0 * J_ * nBond) #J=1.0 tested only

    spin_prop = copy(spin_left)  #the propagated spin state
    for (n, op) in enumerate(operator_list)  #size of the operator list

        if issiteoperator(op) && !isdiagonal(op)
            spin_prop[op.i] ⊻= 1 #spinflip
        else
            while true
                rr = rand()
                if rr < P_h #probability to choose a single-site operator
                    operator_list[n] = SiteOperator{Diagonal}(rand(1:nSpin), lattice)
                    break
                else
                    bond = rand(1:nBond)
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
    for (n, op) in enumerate(operator_list)  #size of the operator list

        if issiteoperator(op)
            site = op.i
            #lower or left leg
            push!(LinkList, First[site])
            push!(LegType, spin_prop[site]) #the spin of the vertex leg
            current_link = length(LinkList)

            if !isdiagonal(op)  #off-diagonal site operator
                spin_prop[site] ⊻= 1 #spinflip
            end

            LinkList[First[site]] = current_link #completes backwards link
            First[site] = current_link + 1
            push!(Associates, nullt)

            #upper or right leg
            push!(LinkList, -99) #we don't yet know what this links to
            push!(LegType, spin_prop[site]) #the spin of the vertex leg
            push!(Associates, nullt)

        else  #diagonal bond operator
            #lower left
            site1 = op.i
            push!(LinkList, First[site1])
            push!(LegType, spin_prop[site1]) #the spin of the vertex leg
            current_link = length(LinkList)
            LinkList[First[site1]] = current_link #completes backwards link
            First[site1] = current_link + 2
            vertex1 = current_link
            push!(Associates, (vertex1 + 1, vertex1 + 2, vertex1 + 3))
            #lower right
            site2 = op.j
            push!(LinkList, First[site2])
            push!(LegType, spin_prop[site2]) #the spin of the vertex leg
            current_link = length(LinkList)
            LinkList[First[site2]] = current_link #completes backwards link
            First[site2] = current_link + 2
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
        current_link = length(LinkList)
        LinkList[First[i]] = current_link
        push!(Associates, nullt)
    end #i

    #DEBUG
    if spin_prop != spin_right
        @debug "Basis state propagation error: LINKED LIST"
    end

end #LinkedList

#############################################################################

#ClusterUpdate
function ClusterUpdate(operator_list, lattice, nSpin, spin_left, spin_right, Associates, LinkList, LegType)

    lsize = length(LinkList)

    #lsize is the size of the linked list
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

    #DEBUG
    #for i = 1:lsize[1]
    #    println(i," ",LegType[i]," ",in_cluster[i])
    #end

    #map back basis states and operator list
    for i in 1:nSpin
        spin_left[i] = LegType[i]  #left basis state
    end

    ocount = nSpin + 1  #next on is leg nSpin + 1
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

    for i in 1:nSpin
        spin_right[i] = LegType[lsize-nSpin+i]  #left basis state
    end

end #ClusterUpdate
