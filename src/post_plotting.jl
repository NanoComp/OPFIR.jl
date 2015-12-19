function post_plotting(p, sol, layer)
    subplot(2,2,1)
    title("non-thermal population for rotational levels")
    hold
    index_offset = p.layer_unknown*(layer-1)
    plot(p.f_dist_ctr, sol[(1:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    plot(p.f_dist_ctr, sol[(2:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    plot(p.f_dist_ctr, sol[(3:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    plot(p.f_dist_ctr, sol[(11:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    plot(p.f_dist_ctr, sol[(12:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    plot(p.f_dist_ctr, sol[(13:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    legend(["L-1", "L", "L+1", "U-1", "U", "U+1"])

    subplot(2,2,2)
    ng0_J4 = (p.ntotal*p.f_G_0/2) * p.C4L .* p.gauss_dist
    ng_J4 = (sol[p.n_rot*p.num_freq + 1 + index_offset] + p.ntotal*p.f_G_0/2) *
            p.C4L .* p.gauss_dist
    plot(p.f_dist_ctr, ng0_J4)
    hold
    plot(p.f_dist_ctr, sol[(2:p.n_rot:p.n_rot*p.num_freq) + index_offset] + ng_J4)
    xlabel("frequency, or velocity")
    ylabel("population (m^-3)")
    legend(["no pump", "pump = $(p.power) W"])
    title("Population of Ground state J=4, K=3")

    subplot(2,2,3)
    n30_J5 = p.ntotal*p.f_3_0/2 * p.C5U .* p.gauss_dist
    n3_J5 = (p.ntotal*p.f_3_0/2 + sol[p.n_rot*p.num_freq + 2 + index_offset]) *
            p.C5U .* p.gauss_dist
    n30_J4 = p.ntotal*p.f_3_0/2 * p.C4U .* p.gauss_dist
    n3_J4 = (p.ntotal*p.f_3_0/2 + sol[p.n_rot*p.num_freq + 2 + index_offset]) *
            p.C4U .* p.gauss_dist
    plot(p.f_dist_ctr, n30_J5)
    hold
    plot(p.f_dist_ctr, sol[(12:p.n_rot:p.n_rot*p.num_freq) + index_offset] + n3_J5) #
    plot(p.f_dist_ctr, (sol[(11:p.n_rot:p.n_rot*p.num_freq)+ index_offset] + n3_J4)*p.g_U/p.g_L) #
    xlabel("frequency, or velocity")
    ylabel("population (m^-3)")
    legend(["no pump", "pump = $(p.power) W"])
    title("Population of V3 state J=5, K=3")

    subplot(2,2,4)

    if p.model_flag==1
        semilogy(p.ntotal/2*[p.f_G_0, p.f_3_0, p.f_6_0, p.f_23, p.f_36, p.f_26,
                             p.f_G_0, p.f_3_0, p.f_6_0, p.f_23, p.f_36, p.f_26],"o-")
        hold
        semilogy(sol[index_offset+p.n_rot*p.num_freq+1:p.n_rot*p.num_freq+p.n_vib+index_offset] +
                 p.ntotal/2*[p.f_G_0, p.f_3_0, p.f_6_0, p.f_23, p.f_36, p.f_26,
                             p.f_G_0, p.f_3_0, p.f_6_0, p.f_23, p.f_36, p.f_26],"o-")

    else
        semilogy(p.ntotal/2*[p.f_G_0, p.f_3_0, p.f_6_0,
                             p.f_G_0, p.f_3_0, p.f_6_0],"o-")
        hold
        semilogy(sol[index_offset+p.n_rot*p.num_freq+1:p.n_rot*p.num_freq+p.n_vib+index_offset] +
                 p.ntotal/2*[p.f_G_0, p.f_3_0, p.f_6_0,
                             p.f_G_0, p.f_3_0, p.f_6_0],"o-")
    end
    xlabel("thermal pools")

    return sol[index_offset + p.n_rot*p.num_freq+1:p.layer_unknown + index_offset]
            #+ p.ntotal/2*[p.f_G_0, p.f_3_0, p.f_6_0, p.f_23, p.f_36, p.f_26]
end

# function plot_thermalpool_layer(p, sol, layer)
#     ######################## plotting and postprocessing ###########
#
#     layer_unknown::Int64 = p.n_rot*p.num_freq + n_vib
#     index_offset = layer_unknown*(layer-1)
#     semilogy(p.ntotal/2*[p.f_G_0, p.f_3_0, p.f_6_0, p.f_23, p.f_36, p.f_26, p.f_G_0, p.f_3_0, p.f_6_0],"o-")
#     hold
#     semilogy(sol[index_offset+p.n_rot*p.num_freq+1:p.n_rot*p.num_freq+9+index_offset] +
#     p.ntotal/2*[p.f_G_0, p.f_3_0, p.f_6_0, p.f_23, p.f_36, p.f_26, p.f_G_0, p.f_3_0, p.f_6_0],"o-")
#     xlabel("thermal pools")
#     # legend(["no pump", "pump = $power W"])
#
#     return sol[index_offset + p.n_rot*p.num_freq+1:layer_unknown + index_offset] #+
#     # p.ntotal/2*[p.f_G_0, p.f_3_0, p.f_6_0, p.f_23, p.f_36, p.f_26, p.f_G_0, p.f_3_0]
# end


# function plot_L_vs_freq_layer(p::pmeter, sol::Array, pressure, power, layer::Int64)
#         layer_unknown::Int64 = p.n_rot*p.num_freq + p.n_vib
#         index_offset = layer_unknown*(layer-1)
#
#         ng0_J4 = (p.ntotal*p.f_G_0/2) * p.C4L .* p.gauss_dist
#         ng_J4 = (sol[18*p.num_freq + 1 + index_offset] + p.ntotal*p.f_G_0/2) *
#                 p.C4L .* p.gauss_dist
#         plot(p.f_dist_ctr, ng0_J4)
#         hold
#         plot(p.f_dist_ctr, sol[(2:p.n_rot:p.n_rot*p.num_freq) + index_offset] + ng_J4)
#         xlabel("frequency, or velocity")
#         ylabel("population (m^-3)")
#         legend(["no pump", "pump = $power W"])
#         title("Population of Ground state J=4, K=3")
#
#         return
# end
#
# function plot_U_vs_freq_layer(p::pmeter, sol::Array, pressure, power, layer::Int64)
#     layer_unknown::Int64 = p.n_rot*p.num_freq + n_vib
#     index_offset = layer_unknown*(layer-1)
#
#     n30_J5 = p.ntotal*p.f_3_0/2 * p.C5U .* p.gauss_dist
#     n3_J5 = (p.ntotal*p.f_3_0/2 + sol[p.n_rot*p.num_freq + 2 + index_offset]) *
#             p.C5U .* p.gauss_dist
#     n30_J4 = p.ntotal*p.f_3_0/2 * p.C4U .* p.gauss_dist
#     n3_J4 = (p.ntotal*p.f_3_0/2 + sol[p.n_rot*p.num_freq + 2 + index_offset]) *
#             p.C4U .* p.gauss_dist
#     plot(p.f_dist_ctr, n30_J5)
#     hold
#     plot(p.f_dist_ctr, sol[(12:p.n_rot:p.n_rot*p.num_freq) + index_offset] + n3_J5)
#     # plot(p.f_dist_ctr, sol[(11:p.n_rot:p.n_rot*p.num_freq) + index_offset] + n3_J4)
#     xlabel("frequency, or velocity")
#     ylabel("population (m^-3)")
#     legend(["no pump", "pump = $power W"])
#     title("Population of V3 state J=5, K=3")
#
#         return
# end
#
# function plot_inv_U(p::pmeter, sol::Array, layer::Int64)
#     pressure = p.pressure
#     freq = p.f_dist_ctr
#     if kvs==0
#         title("U inv, w/o bi collision, pressure = $pressure")
#     else
#         title("U inv, w bi collision, pressure = $pressure")
#     end
#     plot(freq, get_inv_U(p, sol, layer))
#     # legend(["layer number = $layer"])
# end
#
# function plot_inv_L(p::pmeter, sol::Array, layer::Int64)
#     pressure = p.pressure
#     freq = p.f_dist_ctr
#     if kvs==0
#         title("L inv, w/o bi collision, pressure = $pressure")
#     else
#         title("L inv, w bi collision, pressure = $pressure")
#     end
#     plot(freq, get_inv_L(p, sol, layer))
#     # legend(["layer number = $layer"])
# end
#
# function plot_f0()
#     plot([f₀, f₀], [0, 1e15], "r--")
#     hold
#     plot([f_pump, f_pump], [0, 1e15], "r--")
# end
