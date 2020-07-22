# # TODO: This method should be generalizes and moved to the back-end

function copy_stack_down!(
    dg::DGModel,
    m::AtmosModel{FT},
    turbconv::EDMF{FT},
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
) where {FT}
    kernel_calls[:copy_stack_down!] = true
    N_up = n_updrafts(m.turbconv)
    vs_a_i = vars_state(m, Auxiliary(), UInt16)
    vim_a = Vars{vs_a_i}(1:varsize(vs_a_i))

    for i_up in 1:N_up
        i_H = vim_a.turbconv.updraft[i_up].H
        aux[:,i_H, :] .= aux[end,i_H, end]
    end

end
