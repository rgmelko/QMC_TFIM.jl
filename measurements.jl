# measurements.jl
#
# Defines estimators and provides measurements

function Measure()

    spin_prop = copy(spin_left)  

    for i = 1:M #propagate half the list only (to the middle)  
       if operator_list[i,1] == -2
           spin_prop[operator_list[i,2]] = xor(spin_prop[operator_list[i,2]],1) #spinflip
       end
    end #for

    mag = 0
    for j=1:nSpin
        mag += 2*spin_prop[j] - 1
    end

	return mag*mag

end #M2

