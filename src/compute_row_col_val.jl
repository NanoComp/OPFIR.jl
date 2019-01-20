function put_row_col_val(p, rowind, colind, value, row, col, val, s)
    lastrows = vcat(p.layer_unknown:p.layer_unknown:p.layer_unknown*p.num_layers)
    if row in lastrows
        return s
    else
        push!(rowind, row)
        push!(colind, col)
        push!(value, val)
        return (s+1)
    end
        # push!(rowind, row)
        # push!(colind, col)
        # push!(value, val)
        # return (s+1)
end

function compute_row_col_val(rowind, colind, value, p, sol_0)
    s = 1
    # rotational levels in V0 and V3
    tic()
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
            for li in 1:p.n_rot
                row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + li
                for lj in vcat(li-1, li+1) #1:p.n_rot
                    if lj>0 && lj<p.n_rot+1 && p.kDDmat[li, lj] > 0
                        col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + lj
                        ## D-D collision
                        s = put_row_col_val(p, rowind, colind, value, row, row, -p.kDDmat[li,lj], s)
                        s = put_row_col_val(p, rowind, colind, value, row, col, p.kDDmat[lj,li], s)
                    end
                end
                ## SPT term
                s = put_row_col_val(p, rowind, colind, value, row, row, -p.ka[1], s)
            end
        end
    end
    toc()
    ### pump term:
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
            ### pumping L:
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+1)
            val = - p.pump_IR[vi, ri]
            s = put_row_col_val(p, rowind, colind, value, row, row, val, s)

            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
            val = p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(p, rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
            val = - p.pump_IR[vi, ri] * p.CL * p.gauss_dist[vi]
            s = put_row_col_val(p, rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 7
            val = p.pump_IR[vi, ri] * p.g_L/p.g_U * p.CU * p.gauss_dist[vi]
            s = put_row_col_val(p, rowind, colind, value, row, col, val, s)

            ### pumping for U level
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+1) #
            val = p.pump_IR[vi, ri]
            s = put_row_col_val(p, rowind, colind, value, row, col, val, s)

            col = row
            val = - p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(p, rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
            val = p.pump_IR[vi, ri] * p.CL * p.gauss_dist[vi]
            s = put_row_col_val(p, rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 7
            val = - p.pump_IR[vi, ri] * p.g_L/p.g_U * p.CU * p.gauss_dist[vi]
            s = put_row_col_val(p, rowind, colind, value, row, col, val, s)
        end
    end

    ### stimulated emission term:
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+1) #V0, J4 level
            val = - p.WiL*(p.g_L+2)/p.g_L ## J4 -> J5 (diagonal):
            s = put_row_col_val(p, rowind, colind, value, row, row, val, s)
            val = p.WiL ## J5 -> J4:
            s = put_row_col_val(p, rowind, colind, value, row, row+1, val, s)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+2) #V0, J5 level
            val = - p.WiL  # J5 -> J4 (diagonal)
            s = put_row_col_val(p, rowind, colind, value, row, row, val, s)
            val = p.WiL*(p.g_L+2)/p.g_L # J4 -> J5
            s = put_row_col_val(p, rowind, colind, value, row, row-1, val, s)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0) #V3, J4
            val = - p.WiU*p.g_U/(p.g_U-2) ## J4 -> J5 (diagonal)
            s = put_row_col_val(p, rowind, colind, value, row, row, val, s)
            val = p.WiU ## J5 -> J4:
            s = put_row_col_val(p, rowind, colind, value, row, row+1, val, s)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1) #V3, J5
            val = - p.WiU ## J5 -> J4
            s = put_row_col_val(p, rowind, colind, value, row, row, val, s)
            val = p.WiU*p.g_U/(p.g_U-2) ## J4 -> J5
            s = put_row_col_val(p, rowind, colind, value, row, row-1, val, s)
        end
    end
    
    ### transitions between Np
    for ri in 1:p.num_layers
        for pi in 1:p.n_vib
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + pi
            for pj in 1:p.n_vib
                if p.kVVmat[pi, pj] > 1e-8
                    col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + pj
                    s = put_row_col_val(p, rowind, colind, value, row, row, -p.kVVmat[pi, pj], s)
                    s = put_row_col_val(p, rowind, colind, value, row, col, p.kVVmat[pj, pi], s)
                end
            end
        end
    end
    ### SPT from rot to vib levels
    for ri in 1:p.num_layers
        # V0:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
        for vi in 1:p.num_freq
            for k in 1:p.n_rot÷2
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + k
                val = p.ka[1]
                s = put_row_col_val(p, rowind, colind, value, row, col, val, s)
            end
        end
        # V3:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 7
        for vi in 1:p.num_freq
            for k in 1:p.n_rot÷2
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot+ p.n_rot÷2 + k
                val = p.ka[1]
                s = put_row_col_val(p, rowind, colind, value, row, col, val, s)
            end
        end
    end
    ### diffusion for rot/vib levels
    for ri in 2:p.num_layers-1
        index_diffu = vcat(1:p.layer_unknown) + p.layer_unknown * (ri-1)
        for k in index_diffu
            row = k
            col = row - p.layer_unknown
            val = p.D*(-1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
            s = put_row_col_val(p, rowind, colind, value, row, col, val, s)

            val = p.D*(-2.0/(p.Δr)^2)
            s = put_row_col_val(p, rowind, colind, value, row, row, val, s)

            val = p.D*(1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
            col = row + p.layer_unknown
            s = put_row_col_val(p, rowind, colind, value, row, col, val, s)
        end
    end
    # first layer, Neumann BC, for both rotational and vibrational levels
    index_diffu = vcat(1:p.layer_unknown)
    for k in index_diffu
        row = k
        val = p.D*(-1.0/(p.Δr)^2 - 1.0/2.0/p.Δr/p.r_int[1])
        s = put_row_col_val(p, rowind, colind, value, row, row, val, s)

        col = row + p.layer_unknown
        val = p.D*(1.0/(p.Δr)^2 + 1.0/2/p.Δr/p.r_int[1])
        s = put_row_col_val(p, rowind, colind, value, row, col, val, s)
    end
    # BC for the wall: Neumann BC for rotational levels
    Ri = p.num_layers
    index_diffu = p.layer_unknown * (Ri-1) + 1 : p.layer_unknown * Ri - p.n_vib
    for k in index_diffu
        row = k
        col = row - p.layer_unknown
        val = p.D*(-1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
        s = put_row_col_val(p, rowind, colind, value, row, col, val, s)

        val = p.D*(-2.0/(p.Δr)^2) + p.D*(1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
        s = put_row_col_val(p, rowind, colind, value, row, row, val, s)
    end


    # vib level for layer Ri:
    for k in p.layer_unknown*Ri-p.n_vib+1:p.layer_unknown*Ri
        row = k
        col = row - p.layer_unknown
        val = p.D*(-1.0/(2*p.Δr*p.r_int[Ri]) + 1.0/(p.Δr)^2)
        s = put_row_col_val(p, rowind, colind, value, row, col, val, s)

        val = p.D*(-2.0/(p.Δr)^2)
        s = put_row_col_val(p, rowind, colind, value, row, row, val, s)

        val = p.D*(1.0/(2*p.Δr*p.r_int[Ri]) + 1.0/(p.Δr)^2)
        col = row + p.n_vib
        s = put_row_col_val(p, rowind, colind, value, row, col, val, s)
    end
        
    # vib level for layer Ri+1, Robin BC:
    vbar = p.v_avg/2/sqrt(2)/p.norm_time
    for k in 1:p.n_vib
        row = p.layer_unknown*Ri + k
        val = p.D/p.Δr + vbar/4
        s = put_row_col_val(p, rowind, colind, value, row, row, val, s)
        val = -p.D/p.Δr + vbar/4
        s = put_row_col_val(p, rowind, colind, value, row, row-p.n_vib, val, s)
        for pprime in 1:p.n_vib
            val = -vbar * p.frel_vib[k]/4
            col = p.layer_unknown*Ri + pprime
            s = put_row_col_val(p, rowind, colind, value, row, col, val, s)
            s = put_row_col_val(p, rowind, colind, value, row, col-p.n_vib, val, s)
        end
    end
end

# function totN0r(p, sol, ri)
#     tot = 0.
#     for vi in 1:p.num_freq
#         rowoffset = p.layer_unknown*(ri-1) + (vi-1)*p.n_rot
#         tot += sum(sol[rowoffset+1:rowoffset+p.n_rot÷2])
#     end
#     return tot + p.ntotal*p.f_G_0
# end

# function totN3r(p, sol, ri)
#     tot = 0.
#     for vi in 1:p.num_freq
#         rowoffset = p.layer_unknown*(ri-1) + (vi-1)*p.n_rot + p.n_rot÷2
#         tot += sum(sol[rowoffset+1:rowoffset+p.n_rot÷2])
#     end
#     return tot + p.ntotal*p.f_3_0
# end
