function compute_rhs(rhs, p)

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
end
