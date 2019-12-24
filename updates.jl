# updates.jl
#
# Defines the functions that build the linked list and operator cluster update

############################ FUNCTIONS ######################################


lsize = 0
nullt = (0,0,0) #a null tuple

#LinkedList
function LinkedList()

    #initialize linked list data structures
    LinkList = zeros(Int,0)  #needed for cluster update
    LegType = zeros(Int,0)  

    #A diagonal bond operator has non trivial associates for cluster building
    #Associates = zeros((Int,Int,Int),0)
    global Associates = []

    First = zeros(Int,0)  #scope is this function only

    #The first N elements of the linked list are the spins of the LHS basis state
    for i = 1:nSpin
        push!(First,i)
        push!(LinkList,-99) #we don't yet know what these link to
        push!(LegType,spin_left[i])
        push!(Associates,nullt) #no nontrivial associates for train wf spins
    end #i
    println(Associates)

    spin_prop = copy(spin_left)  #the propagated spin state

    #Now, add the 2M operators to the linked list.  Each has either 2 or 4 legs
    for i = 1:2*M  #size of the operator list

       if operator_list[i,1] == -2  #off-diagonal site operator
           site = operator_list[i,2]
           spin_prop[i] = xor(spin_prop[i],1) #spinflip

       elseif operator_list[i,1] == -1  #diagonal site operator
           site = operator_list[i,2]
           #lower or left leg
           push!(LinkList,First[site])
           push!(LegType,spin_prop[site]) #the spin of the vertex leg
           current_link = size(LinkList)
           LinkList[First[site]] = current_link[1] #completes backwards link
           First[site] = current_link[1] + 1
           push!(Associates,nullt)
           #upper or right leg
           push!(LinkList,-99) #we don't yet know what this links to
           push!(LegType,spin_prop[site]) #the spin of the vertex leg
           push!(Associates,nullt)

       else  #diagonal bond operator
           #lower left
           site1 = operator_list[i,1]
           push!(LinkList,First[site1])
           push!(LegType,spin_prop[site1]) #the spin of the vertex leg
           current_link = size(LinkList)
           LinkList[First[site1]] = current_link[1] #completes backwards link
           First[site1] = current_link[1] + 2
           vertex1 = current_link[1]
           push!(Associates,(vertex1+1,vertex1+2,vertex1+3))
           #lower right
           site2 = operator_list[i,2]
           push!(LinkList,First[site2])
           push!(LegType,spin_prop[site2]) #the spin of the vertex leg
           current_link = size(LinkList)
           LinkList[First[site2]] = current_link[1] #completes backwards link
           First[site2] = current_link[1] + 2
           push!(Associates,(vertex1,vertex1+2,vertex1+3))
           #upper left
           push!(LinkList,-99) #we don't yet know what this links to
           push!(LegType,spin_prop[site1]) #the spin of the vertex leg
           push!(Associates,(vertex1,vertex1+1,vertex1+3))
           #upper right 
           push!(LinkList,-99) #we don't yet know what this links to
           push!(LegType,spin_prop[site2]) #the spin of the vertex leg
           push!(Associates,(vertex1,vertex1+1,vertex1+2))

       end #if

    end #i

    #The last N elements of the linked list are the final spin state
    for i = 1:nSpin
        push!(LinkList,First[i]) 
        push!(LegType,spin_prop[i]) 
        current_link = size(LinkList)
        LinkList[First[i]] = current_link[1]
        push!(Associates,nullt)
    end #i

    global lsize = size(LinkList)

    #DEBUG
    println(lsize)
    for i = 1:lsize[1]
       println(i," ",LinkList[i])
       #println(Associates[i])
       #if Associates[i] == nullt
       #    println("NULL")
       #end
    end
    println(LegType," ",size(LegType))

    #DEBUG
    if spin_prop != spin_right
        println("Basis state propagation error: LINKED LIST")
    end

end #LinkedList

#############################################################################

#ClusterUpdate
function ClusterUpdate()

    #lsize is the size of the linked list
    in_cluster=zeros(Int,lsize[1])

    cstack = zeros(Int,0)  #This is the stack of vertices in a cluster

    for i = 1:lsize[1]

        #Add a new leg onto the cluster
        if (in_cluster[i] == 0 && Associates[i] == nullt)

            println("add new ",i)

        end #if


    end #for i

end #ClusterUpdate


