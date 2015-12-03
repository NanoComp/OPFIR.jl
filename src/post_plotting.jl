function post_plotting(para::parameter, nend::Array, pressure, power, layer::Int64)
    subplot(2,2,1)
    ######################## plotting and postprocessing ###########

    layer_unknown::Int64 = n_rot*para.num_freq + n_vib;

    title("non-thermal population for rotational levels");
    hold;
    index_offset = layer_unknown*(layer-1);
    plot(para.f_dist_ctr, nend[(1:n_rot:n_rot*para.num_freq) + index_offset],"-")
    plot(para.f_dist_ctr, nend[(2:n_rot:n_rot*para.num_freq) + index_offset],"-")
    plot(para.f_dist_ctr, nend[(3:n_rot:n_rot*para.num_freq) + index_offset],"-")
    plot(para.f_dist_ctr, nend[(11:n_rot:n_rot*para.num_freq) + index_offset],"-")
    plot(para.f_dist_ctr, nend[(12:n_rot:n_rot*para.num_freq) + index_offset],"-")
    plot(para.f_dist_ctr, nend[(13:n_rot:n_rot*para.num_freq) + index_offset],"-")
    legend(["L-1", "L", "L+1", "U-1", "U", "U+1"])

    subplot(2,2,2)
    ng0_J4 = (para.ntotal*f_G/2) * C4L .* para.gauss_dist;
    ng_J4 = (nend[18*para.num_freq + 1 + index_offset] + para.ntotal*f_G/2) *
            C4L .* para.gauss_dist;
    plot(para.f_dist_ctr, ng0_J4);
    hold;
    plot(para.f_dist_ctr, nend[(2:n_rot:n_rot*para.num_freq) + index_offset] + ng_J4);
    xlabel("frequency, or velocity");
    ylabel("population (m^-3)")
    legend(["no pump", "pump = $power W"])
    title("Population of Ground state J=4, K=3")

    subplot(2,2,3)
    n30_J5 = para.ntotal*f_3/2 * C5U .* para.gauss_dist;
    n3_J5 = (para.ntotal*f_3/2 + nend[n_rot*para.num_freq + 2 + index_offset]) *
            C5U .* para.gauss_dist;
    n30_J4 = para.ntotal*f_3/2 * C4U .* para.gauss_dist;
    n3_J4 = (para.ntotal*f_3/2 + nend[n_rot*para.num_freq + 2 + index_offset]) *
            C4U .* para.gauss_dist;
    plot(para.f_dist_ctr, n30_J5);
    hold;
    plot(para.f_dist_ctr, nend[(12:n_rot:n_rot*para.num_freq) + index_offset] + n3_J5); #
    plot(para.f_dist_ctr, (nend[(11:n_rot:n_rot*para.num_freq)+ index_offset] + n3_J4)*g_U/g_L) #
    xlabel("frequency, or velocity");
    ylabel("population (m^-3)")
    legend(["no pump", "pump = $power W"])
    title("Population of V3 state J=5, K=3")

    subplot(2,2,4)
#         plot(para.nend[1:18],"o-")
#         xlabel("rotational level label")
    semilogy(para.ntotal/2*[f_G, f_3, f_6, f_23, f_36, f_26, f_G, f_3],"o-")
    hold
    semilogy(nend[index_offset+n_rot*para.num_freq+1:n_rot*para.num_freq+n_vib+index_offset] +
    para.ntotal/2*[f_G, f_3, f_6, f_23, f_36, f_26, f_G, f_3, f_6, f_23, f_36, f_26],"o-")
    xlabel("thermal pools")
#     legend(["no pump", "pump = $power W"])

    return nend[index_offset + n_rot*para.num_freq+1:layer_unknown + index_offset]
            #+ para.ntotal/2*[f_G, f_3, f_6, f_23, f_36, f_26]
end

function plot_thermalpool_layer(para::parameter, nend::Array, pressure, power, layer::Int64)
    ######################## plotting and postprocessing ###########

    layer_unknown::Int64 = n_rot*para.num_freq + n_vib;
    index_offset = layer_unknown*(layer-1);
    semilogy(para.ntotal/2*[f_G, f_3, f_6, f_23, f_36, f_26, f_G, f_3, f_6],"o-")
    hold
    semilogy(nend[index_offset+n_rot*para.num_freq+1:n_rot*para.num_freq+9+index_offset] +
    para.ntotal/2*[f_G, f_3, f_6, f_23, f_36, f_26, f_G, f_3, f_6],"o-")
    xlabel("thermal pools")
    # legend(["no pump", "pump = $power W"])

    return nend[index_offset + n_rot*para.num_freq+1:layer_unknown + index_offset] #+
    # para.ntotal/2*[f_G, f_3, f_6, f_23, f_36, f_26, f_G, f_3]
end


function plot_L_vs_freq_layer(para::parameter, nend::Array, pressure, power, layer::Int64)
        layer_unknown::Int64 = n_rot*para.num_freq + n_vib;
        index_offset = layer_unknown*(layer-1);

        ng0_J4 = (para.ntotal*f_G/2) * C4L .* para.gauss_dist;
        ng_J4 = (nend[18*para.num_freq + 1 + index_offset] + para.ntotal*f_G/2) *
                C4L .* para.gauss_dist;
        plot(para.f_dist_ctr, ng0_J4);
        hold;
        plot(para.f_dist_ctr, nend[(2:n_rot:n_rot*para.num_freq) + index_offset] + ng_J4);
        xlabel("frequency, or velocity");
        ylabel("population (m^-3)")
        legend(["no pump", "pump = $power W"])
        title("Population of Ground state J=4, K=3")

        return;
end

function plot_U_vs_freq_layer(para::parameter, nend::Array, pressure, power, layer::Int64)
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib;
    index_offset = layer_unknown*(layer-1);

    n30_J5 = para.ntotal*f_3/2 * C5U .* para.gauss_dist;
    n3_J5 = (para.ntotal*f_3/2 + nend[n_rot*para.num_freq + 2 + index_offset]) *
            C5U .* para.gauss_dist;
    n30_J4 = para.ntotal*f_3/2 * C4U .* para.gauss_dist;
    n3_J4 = (para.ntotal*f_3/2 + nend[n_rot*para.num_freq + 2 + index_offset]) *
            C4U .* para.gauss_dist;
    plot(para.f_dist_ctr, n30_J5);
    hold;
    plot(para.f_dist_ctr, nend[(12:n_rot:n_rot*para.num_freq) + index_offset] + n3_J5);
    # plot(para.f_dist_ctr, nend[(11:n_rot:n_rot*para.num_freq) + index_offset] + n3_J4);
    xlabel("frequency, or velocity");
    ylabel("population (m^-3)")
    legend(["no pump", "pump = $power W"])
    title("Population of V3 state J=5, K=3")

        return;
end

function plot_inv_U(para::parameter, sol::Array, layer::Int64)
    pressure = para.pressure
    freq = para.f_dist_ctr
    if kvs==0
        title("U inv, w/o bi collision, pressure = $pressure")
    else
        title("U inv, w bi collision, pressure = $pressure")
    end
    plot(freq, get_inv_U(para, sol, layer))
    # legend(["layer number = $layer"])
end

function plot_inv_L(para::parameter, sol::Array, layer::Int64)
    pressure = para.pressure
    freq = para.f_dist_ctr
    if kvs==0
        title("L inv, w/o bi collision, pressure = $pressure")
    else
        title("L inv, w bi collision, pressure = $pressure")
    end
    plot(freq, get_inv_L(para, sol, layer))
    # legend(["layer number = $layer"])
end

function plot_f0()
    plot([f₀, f₀], [0, 1e15], "r--")
    hold
    plot([f_pump, f_pump], [0, 1e15], "r--")
end
