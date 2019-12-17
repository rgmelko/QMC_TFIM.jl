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
	println(First)
	println(LinkList)
	println(LegType)
	println(Associates)

    spin_prop = copy(spin_left)  #the propagated spin state

    #Now, add the 2M operators to the linked list.  Each has either 2 or 4 legs
    for i = 1:2*M  #size of the operator list

       if operator_list[i,1] == -2  #off-diagonal site operator
		   site = operator_list[i,2]
           spin_prop[i] = xor(spin_prop[i],1) #spinflip

       elseif operator_list[i,1] == -1
           println("diagonal site operator at ",i)
	   else
           println("diagonal bond operator at ",i)
	   end #if

	end #i
end #LinkedList

#############################################################################
