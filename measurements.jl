# measurements.jl
#
# Defines estimators and provides measurements
using Statistics


function sample(spin_left, operator_list)
    spin_prop = copy(spin_left)
    for o in operator_list[1:M] # propagate half the list only (to the middle)
        if o isa SiteOperator{OffDiagonal}
            spin_prop[o.i] ⊻= 1 # spinflip
        end
    end
    return spin_prop
end

magnetization(spin_prop) = mean(x->2x - 1, spin_prop)

function energy_abs_zero(h, J, spin_prop, operator_list)
    m_d = count(x->isdiagonal(x) && issiteoperator(x), operator_list)
    # E = (J - h/2)*m_d/M
    # E = J*m_d + h/2*(2M-m_d)
    # E /= M

    # E_J = -2J * m_d / M
    E_J = (2 * M * h / m_d)

    E_h = 0
    inv_n = 1 / m_d
    return E_J, E_h, inv_n, m_d
end
