# updates.jl
#
# Defines the functions that build the linked list and operator cluster update

############################ FUNCTIONS ######################################

lsize = 0
nullt = (0,0,0) #a null tuple

#Diagonal update
function DiagonalUpdate()
    spin_prop = copy(spin_left)  #the propagated spin state
    for i = 1:2*M  #size of the operator list

       if operator_list[i,1] == -2
           spin_prop[operator_list[i,2]] = xor(spin_prop[operator_list[i,2]],1) #spinflip

       else
           flag = false
           while flag == false
               rr = rand() 
               if P_h > rr #probability to choose a single-site operator
                   operator_list[i,1] = -1
                   site = rand(1:nSpin)
                   operator_list[i,2] = site
                   flag = true
                   #println(i," site ",site)
               else
                   bond = rand(1:nBond)
                   if spin_prop[bond_spin[bond,1]] == spin_prop[bond_spin[bond,2]] #spins must be the same
                       operator_list[i,1] = bond_spin[bond,1]
                       operator_list[i,2] = bond_spin[bond,2]
                       flag = true
                       #println(i," bond ",bond)
                   end#if
                   #println(P_h," ",rr," bond")
               end
           end #while

        end#if
    end #for

    #DEBUG
    if spin_prop != spin_right  #check the spin propagation for error
        println("Basis state propagation error: DiagonalUpdate")
    end


end #DiagonalUpdate

#############################################################################

#LinkedList
function LinkedList()

    #initialize linked list data structures
    global LinkList = zeros(Int,0)  #needed for cluster update
    global LegType = zeros(Int,0)  

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

    spin_prop = copy(spin_left)  #the propagated spin state

    #Now, add the 2M operators to the linked list.  Each has either 2 or 4 legs
    for i = 1:2*M  #size of the operator list

       if operator_list[i,1] == -2  #off-diagonal site operator
           site = operator_list[i,2]
           #lower or left leg
           push!(LinkList,First[site])
           push!(LegType,spin_prop[site]) #the spin of the vertex leg
           spin_prop[site] = xor(spin_prop[site],1) #spinflip
           current_link = size(LinkList)
           LinkList[First[site]] = current_link[1] #completes backwards link
           First[site] = current_link[1] + 1
           push!(Associates,nullt)
           #upper or right leg
           push!(LinkList,-99) #we don't yet know what this links to
           push!(LegType,spin_prop[site]) #the spin of the vertex leg
           push!(Associates,nullt)

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
    #println(lsize)
    #for i = 1:lsize[1]
    #   println(i," ",LinkList[i])
    #   #println(Associates[i])
    #   #if Associates[i] == nullt
    #   #    println("NULL")
    #   #end
    #end
    #println(LegType," ",size(LegType))

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

    ccount = 0 #cluster number counter
    for i = 1:lsize[1]

        #Add a new leg onto the cluster
        if (in_cluster[i] == 0 && Associates[i] == nullt)

            ccount+=1
            push!(cstack,i) 
            in_cluster[cstack[end]] = ccount  

            flip = rand(Bool) #flip a coin for the SW cluster flip
            if flip == true 
                LegType[cstack[end]] =  xor(LegType[cstack[end]],1) #spinflip
            end

            while isempty(cstack) == false

                leg = LinkList[cstack[end]]
                pop!(cstack)
                #println("leg ",leg," ",cstack)

                if in_cluster[leg] == 0

                    in_cluster[leg] = ccount; #add the new leg and flip it 
                    if flip == true 
                        #println("gonna flip ",leg)
                        LegType[leg] =  xor(LegType[leg],1) 
                    end
                    #now check all associates and add to cluster
                    assoc = Associates[leg] #a 3-element array
                    if assoc != nullt
                        push!(cstack,assoc[1])
                        in_cluster[assoc[1]] = ccount
                        push!(cstack,assoc[2])
                        in_cluster[assoc[2]] = ccount
                        push!(cstack,assoc[3])
                        in_cluster[assoc[3]] = ccount
                        if flip == true 
                            LegType[assoc[1]] =  xor(LegType[assoc[1]],1) 
                            LegType[assoc[2]] =  xor(LegType[assoc[2]],1) 
                            LegType[assoc[3]] =  xor(LegType[assoc[3]],1) 
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
    for i = 1:nSpin
        spin_left[i] = LegType[i]  #left basis state
        ocount += 1
    end 

    ocount += 1  #next on is leg nSpin + 1
    for i = 1:2*M  
        if operator_list[i,1] != -2 && operator_list[i,1] != -1
            ocount += 4
        else
            if LegType[ocount] == LegType[ocount + 1]  #diagonal
                operator_list[i,1] = -1
                #println("DCHANGE ",i," ",ocount)
            else
                operator_list[i,1] = -2 #off-diagonal
                #println("OCHANGE ",i," ",ocount)
            end
            ocount += 2
        end
    end

    for i = 1:nSpin
        spin_right[i] = LegType[lsize[1] - nSpin + i]  #left basis state
    end 


end #ClusterUpdate


