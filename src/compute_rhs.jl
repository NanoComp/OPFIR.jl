function compute_rhs(rhs, p, sol)
    for ri in 1:p.num_layers
        for vi in 1:p.num_freq

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+1)
            rhs[row] = p.pump_IR[vi, ri] * p.gauss_dist[vi] * p.ntotal *
                       (p.CL * p.f_G_0/2 - p.CU * p.f_3_0/2 * p.g_L/p.g_U)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
            rhs[row] = - p.pump_IR[vi, ri] * p.gauss_dist[vi] * p.ntotal *
                       (p.CL * p.f_G_0/2 - p.CU * p.f_3_0/2 * p.g_L/p.g_U)
        end
    end
end


function pump_total(p, sol, layer)
    L_tot = Nl_total_dist_layer(p, sol, layer)
    U_tot = Nu_total_dist_layer(p, sol, layer)
    pump_dist = p.pumpR .* (L_tot-p.g_L/p.g_U*U_tot)
    return  sum(pump_dist)
end
