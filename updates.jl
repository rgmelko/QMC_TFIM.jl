# updates.jl
#
# Defines the functions that build the linked list and operator cluster update

############################ FUNCTIONS ######################################

#LinkedList
function LinkedList()
    #spin_prop = copy(spin_left)  #the propagated spin state
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
