function compute_rhs(rhs, p, sol)

    for ri in 1:p.num_layers
        for vi in 1:p.num_freq
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 2
            rhs[row] = - p.pumpR[vi] * p.gauss_dist[vi] * p.ntotal *
                       (p.C4L * p.f_G/2 - p.C5U * p.f_3/2 * p.g_L/p.g_U)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + 3
            rhs[row] = p.pumpR[vi] * p.gauss_dist[vi] * p.ntotal *
                       (p.C4L * p.f_G/2 - p.C5U * p.f_3/2 * p.g_L/p.g_U)
        end
    end

    if p.model_flag==2
        for ri in 1:p.num_layers
            # row of V3A level:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
            tot_pump = pump_total(p, sol, ri)
            # V3A -> VΣA net rate
            rhs[row] = tot_pump * sol[row]/(sol[row]+sol[row+1])
            rhs[row+1] = - tot_pump * sol[row]/(sol[row]+sol[row+1])
        end
    end

end


function pump_total(p, sol, layer)
    L_tot = Nl_total_dist_layer(p, sol, layer)
    U_tot = Nu_total_dist_layer(p, sol, layer)
    pump_dist = p.pumpR .* (L_tot-p.g_L/p.g_U*U_tot)
    return  sum(pump_dist)
end
