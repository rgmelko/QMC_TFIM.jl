# updates.jl
#
# Defines the functions that build the linked list and operator cluster update

############################ FUNCTIONS ######################################

#LinkedList
function LinkedList()

    #linked list data structures
    LinkList = zeros(Int,0)  #needed for cluster update
    LegType = zeros(Int,0)  
    Associates = zeros(Int,0)  

	First = zeros(Int,0)  #scope is this function only

    #The first N elements of the linked list are the spins of the LHS basis state
	for i = 1:nSpin
		push!(First,i)
		push!(LinkList,-99) #we don't yet know what these link to
		push!(LegType,spin_left[i])
		#push!(Associates,spin_left[i])
	end #i

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
		   current_link = size(LinkList)
		   LinkList[First[site]] = current_link[1] #completes backwards link
		   First[site] = current_link[1] + 1
           #upper or right leg
		   push!(LinkList,-99) #we don't yet know what this links to

	   else  #diagonal bond operator
           #lower left
		   site1 = operator_list[i,1]
		   push!(LinkList,First[site1])
		   current_link = size(LinkList)
		   LinkList[First[site1]] = current_link[1] #completes backwards link
		   First[site1] = current_link[1] + 2
           #lower right
		   site2 = operator_list[i,2]
		   push!(LinkList,First[site2])
		   current_link = size(LinkList)
		   LinkList[First[site2]] = current_link[1] #completes backwards link
		   First[site2] = current_link[1] + 2
           #upper left
		   push!(LinkList,-99) #we don't yet know what this links to
           #upper right 
		   push!(LinkList,-99) #we don't yet know what this links to

	   end #if

	end #i

    #The last N elements of the linked list are the final spin state
	for i = 1:nSpin
		push!(LinkList,First[i]) 
		current_link = size(LinkList)
		LinkList[First[i]] = current_link[1]
	end #i

    lsize = size(LinkList)
    for i = 1:lsize[1]
       println(i," ",LinkList[i])
	end

end #LinkedList

#############################################################################
