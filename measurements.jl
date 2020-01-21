# measurements.jl
#
# Defines estimators and provides measurements
using Statistics


function sample(spin_left, operator_list)

    spin_prop = copy(spin_left)

    for i in 1:M #propagate half the list only (to the middle)
        if operator_list[i, 1] == -2
            spin_prop[operator_list[i, 2]] ⊻= 1 #spinflip
        end
    end

    return spin_prop
end


magnetization(spin_prop) = mean(x->2x-1, spin_prop)

#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j + 1)

function energy_abs_zero(h, J, spin_prop, operator_list)
    m_d = count(x->x>0, operator_list[:, 1])


    E = (J - h/2)*m_d/M
    return -E/size(operator_list, 1)
end