function compute_rhs(rhs, p, sol)
    for ri in 1:p.num_layers
        for vi in 1:p.num_freq

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+1)
            rhs[row] += p.pump_IR[vi, ri] * p.gauss_dist[vi] * p.ntotal *
                       (p.CL * p.f_G_0 - p.CU * p.f_3_0 * p.g_L/p.g_U)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
            rhs[row] += - p.pump_IR[vi, ri] * p.gauss_dist[vi] * p.ntotal *
                       (p.CL * p.f_G_0 - p.CU * p.f_3_0 * p.g_L/p.g_U)
        end
    end

    # flux = posflux(p, sol)
    vbar = p.v_avg/2/sqrt(2)
    for vi in 1:p.num_freq
        for j in 1:p.n_rot
            row = p.layer_unknown*p.num_layers + (vi-1)*p.n_rot + j
            rhs[row] += vbar/p.norm_time * (sol[row] + sol[row-p.layer_unknown])/2 * 1/2
            for jprime in 1:p.n_rot
                rowprime = p.layer_unknown*p.num_layers + (vi-1)*p.n_rot + jprime
                rhs[row] += - p.rotpopfr[j] * vbar/p.norm_time * (sol[rowprime] + sol[rowprime-p.layer_unknown])/2 * 1/2
            end
            rhs[row] += -p.D/p.Δr * sol[row-p.layer_unknown]
        end
    end

end

function rotpopfracv(vi, p)
    # p.velocity[vi] < 0
    ctot = ci = 0.0
    for k in 1:p.num_freq
        if p.velocity[k] < 0
            ctot += p.gauss_dist[k]
            if k == vi
                ci = p.gauss_dist[k]
            end
        end
    end
    return ci/ctot
end

function posflux(p, sol)
    flux = 0.0
    for lprime in 1:p.n_rot
        for viprime in 1:p.num_freq
            if p.velocity[viprime] >= 0
                n = p.layer_unknown*p.num_layers + (viprime-1)*p.n_rot + lprime
                flux += p.velocity[viprime] /2 * (sol[n] + sol[n-p.layer_unknown])
                # flux += p.v_avg /2 * (sol[n] + sol[n-p.layer_unknown])
            end
        end
    end
    return flux
end

#
# function pump_total(p, sol, layer)
#     L_tot = Nl_total_dist_layer(p, sol, layer)
#     U_tot = Nu_total_dist_layer(p, sol, layer)
#     pump_dist = p.pumpR .* (L_tot-p.g_L/p.g_U*U_tot)
#     return  sum(pump_dist)
# end
