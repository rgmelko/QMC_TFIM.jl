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
	for i = 1:nSpin
		push!(First,i)
		push!(LinkList,-99)
		push!(LegType,spin_left[i])
		#push!(Associates,spin_left[i])
	end #i
	println(First)
	println(LinkList)
	println(LegType)
	println(Associates)

    spin_prop = copy(spin_left)  #the propagated spin state

    for i = 1:2*M  #size of the operator list

       if operator_list[i,1] == -2
           println("off-diagonal bond operator at ",i)
       elseif operator_list[i,1] == -1
           println("diagonal site operator at ",i)
	   else
           println("diagonal bond operator at ",i)
	   end #if

	end #i
end #LinkedList

#############################################################################
