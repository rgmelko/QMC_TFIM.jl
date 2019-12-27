# measurements.jl
#
# Defines estimators and provides measurements

function Measure()

    #spin_prop = copy(spin_left)  
    #m2 = 0

    #count = 0
    #for i = 1:2*M  

    #   if operator_list[i,1] == -2
    #       spin_prop[operator_list[i,2]] = xor(spin_prop[operator_list[i,2]],1) #spinflip
    #   end

    #   if mod(i,nSpin) == 0  #the basis state should be roughly independent every N propagation steps
    #       count += 1
    #       mag = 0
    #       for j=1:nSpin
    #           mag += 2*spin_prop[j] - 1
    #       end
    #   m2 += mag*mag
    #   end

    #end #for

    mag = 0
    for j=1:nSpin
        mag += 2*spin_left[j] - 1
    end
	return mag*mag

    #m2 /= count
    #return m2 
    #println(count," ",M1)

end #M2

