# measurements.jl
#
# Defines estimators and provides measurements

function sample(spin_left, operator_list)

    spin_prop = copy(spin_left)

    for i in 1:M #propagate half the list only (to the middle)
        if operator_list[i, 1] == -2
            spin_prop[operator_list[i, 2]] ⊻= 1 #spinflip
        end
    end

    return spin_prop
end


magnetization(spin_prop) = mapreduce(x->2*x - 1, +, spin_prop)
magnetization_sqr(spin_prop) = magnetization(spin_prop) ^ 2


#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j + 1)

function energy_abs_zero(h, J, spin_prop, operator_list)
    m_d, m_o = 0., 0.
    n_h, n_J = 0, 0

    for i in 1:div(size(operator_list, 1), 2)
        if operator_list[i, 1] == -2
            m_o += 1
        else
            m_d += 1
            if operator_list[i, 1] == -1
                n_h += 1
            elseif operator_list[i, 1] > 0
                n_J += 1
            end
        end
    end

    E = m_d + m_o/2
    E -= h*n_h
    E -= J*n_J

    return -E
end