function OPFIR_compute_rhs(rhs::Array, para)
    ######################### right hand side: from pumping term #########################
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib;
    row::Int64 = 0;
    for ri::Int64 in 1:para.num_layers
        for vi::Int64 in 1:para.num_freq
            row = (ri-1)*layer_unknown + (vi-1)*n_rot + 2;
            rhs[row] = - para.pumpR[vi] * para.gauss_dist[vi] * para.ntotal *
                       (C4L * f_G/2 - C5U * f_3/2 * g_L/g_U);
            
            row = (ri-1)*layer_unknown + (vi-1)*n_rot + n_rot/2 + 3;
            rhs[row] = para.pumpR[vi] * para.gauss_dist[vi] * para.ntotal *
                       (C4L * f_G/2 - C5U * f_3/2 * g_L/g_U);
        end
    end
end