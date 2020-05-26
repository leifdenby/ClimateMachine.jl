# # TODO: This method should be generalizes and moved to the back-end

# function compute_updraft_top!(
#     dg::DGModel,
#     m::AtmosModel{FT},
#     turbconv::EDMF,
#     Q::MPIStateArray,
#     t::Real,
#     elems::UnitRange,
# ) where {FT, N}

#     grid = dg.grid
#     vars = vars_state_conservative(m, FT)
#     vrange = 1:size(Q, 3),
#     i::Int = 1
#     j::Int = 1

#     Nq = N + 1
#     Nqk = dimensionality(grid) == 2 ? 1 : Nq

#     var_names = flattenednames(vars)
#     var_ind = findfirst(s -> s == var, var_names)

#     state_data = array_device(Q) isa CPU ? Q.realdata : Array(Q.realdata)
#     z = vrange.start
#     result = state_data[1, var_ind, z]
#     for i in vrange
#         for j in vrange
#             for ev in vrange
#                 for k in 1:Nqk
#                     ijk = i + Nq * ((j - 1) + Nq * (k - 1))
#                     for i_up in 1:n_updrafts(m.turbconv)

#                         vars = Vars{st}(Q[ijk, :, ev])
#                         up_i_ρa = vars.turbconv.updraft[i_up].ρa
#                         ρinv = 1/vars.ρ

#                         if up_i_ρa*ρinv>0
#                             aux_vars["turbconv.updraft[$i].updraft_top"][j] = z[j]
#                         end

#                         new_result = op(result, state_data[ijk, var_ind, ev])
#                         if !isequal(new_result, result)
#                             result = new_result
#                             z = ev
#                             Q[ijk, i_var, ev]
#                         end
#                 end
#             end
#         end
#     end

#     vsc = vars_state_conservative(m, FT)
#     for i in 1:N
#         for j in 1:length(z)
#             up_i_ρa = state_vars["turbconv.updraft[$i].ρa"][j]
#             ρinv = 1/state_vars["ρ"][j]
#             if up_i_ρa*ρinv>0
#                 aux_vars["turbconv.updraft[$i].updraft_top"][j] = z[j]
#             end
#         end
#         # i_updraft_top = varsindex(vsc, :updraft_top)
#         # Q[ijk, i_updraft_top, elems] = aux_vars["turbconv.updraft[$i].updraft_top"][j]
#     end
#     #   -------------  Compute upd_top
#     # YAIR - check this with Charlie
#     # results = []
#     # for i in 1:N
#     #     res = reduce_element_stack(
#     #         x -> max(x, 0),
#     #         grid,
#     #         Q,
#     #         vars_state_conservative(m, FT),
#     #         "turbconv.updraft[$i].ρa";
#     #     )
#     #     push!(results, res)
#     # end
#     # for i in 1:N
#     #     # i_updraft_top = varsindex(vsc, :updraft_top)
#     #     Q[ijk, i_updraft_top, elems] = aux_vars["turbconv.updraft[$i].updraft_top"][j]
#     #     for j in 1:length(up[i].ρa)
#     #         if up[i].ρa*ρinv>eps(FT)
#     #             up_a[i].updraft_top = results[i][]
#     #         end
#     #     end
#     # end

# end
