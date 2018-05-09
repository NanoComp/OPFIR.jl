function put_row_col_val(rowind, colind, value, row, col, val, s)
    rowind[s] = row
    colind[s] = col
    value[s] = val
    return (s+1)
end

function compute_row_col_val(rowind, colind, value, p, sol_0)
    s = 1
    # rotational levels in V0
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
        #     for j in 1:p.n_rot
        #         row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + j
        #         val = -p.ka[j] - p.kro # SPT and V swap
        #         s = put_row_col_val(rowind, colind, value, row, row, val, s)
        #     end

            for li in 1:p.n_rot
                for lj in 1:p.n_rot
                    if p.kDDmat[li, lj] > 0
                        row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + li
                        col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + lj
                        ## D-D collision
                        s = put_row_col_val(rowind, colind, value, row, row, -p.kDDmat[li,lj], s)
                        s = put_row_col_val(rowind, colind, value, row, col, p.kDDmat[lj,li], s)
                    end
                end
            end
        end
    end
    #
    ### pump term:
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
            ### pumping L:
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+1)
            val = - p.pump_IR[vi, ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            if p.n_vib > 0
                col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1 #V0A
                val = - p.pump_IR[vi, ri] * p.CL * p.gauss_dist[vi]
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end

            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
            val = p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
            if p.n_vib > 0
                col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2 #V3A
                val = p.pump_IR[vi, ri] * p.CU * p.gauss_dist[vi] * p.g_L/p.g_U
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
            ### pumping for U level
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+1) #
            val = p.pump_IR[vi, ri]
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
            if p.n_vib > 0
                col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1 # V0A
                val = p.pump_IR[vi, ri] * p.CL * p.gauss_dist[vi]
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end

            col = row
            val = - p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
            if p.n_vib > 0
                col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2 # V3A
                val = - p.pump_IR[vi, ri] * p.CU * p.gauss_dist[vi] * p.g_L/p.g_U
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
        end
    end

    ### stimulated emission term:
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+1) #V0, J4 level
            val = - p.WiL*(p.g_L+2)/p.g_L ## J4 -> J5 (diagonal):
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = p.WiL ## J5 -> J4:
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+2) #V0, J5 level
            val = - p.WiL  # J5 -> J4 (diagonal)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = p.WiL*(p.g_L+2)/p.g_L # J4 -> J5
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0) #V3, J4
            val = - p.WiU*p.g_U/(p.g_U-2) ## J4 -> J5 (diagonal)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = p.WiU ## J5 -> J4:
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1) #V3, J5
            val = - p.WiU ## J5 -> J4
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = p.WiU*p.g_U/(p.g_U-2) ## J4 -> J5
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end
    end

    ########################## add the diffusion term ##########################

    for ri in 2:p.num_layers
        index_diffu = vcat(p.layer_unknown*(ri-1)+1 : p.layer_unknown*ri)
        for k in index_diffu
            row = k
            col = row - p.layer_unknown
            val = p.D*(-1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            val = p.D*(-2.0/(p.Δr)^2)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = p.D*(1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
            col = row + p.layer_unknown
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
        end
    end
    # first layer, Neuemann BC
    index_diffu = vcat(1:p.layer_unknown)
    for k in index_diffu
        row = k
        val = p.D*(-1.0/(p.Δr)^2 - 1.0/2.0/p.Δr/p.r_int[1])
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        col = row + p.layer_unknown
        val = p.D*(1.0/(p.Δr)^2 + 1.0/2/p.Δr/p.r_int[1])
        s = put_row_col_val(rowind, colind, value, row, col, val, s)
    end
    #
    # ################ question: BC for the wall? #######
    Ri = p.num_layers
    # for ri in Ri
    #     index_diffu = vcat(p.layer_unknown*(ri-1)+1 : p.layer_unknown*ri)
    #     for k in index_diffu
    #         row = k
    #         col = row - p.layer_unknown
    #         val = p.D*(-1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
    #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
    #
    #         val = p.D*(-2.0/(p.Δr)^2)
    #         s = put_row_col_val(rowind, colind, value, row, row, val, s)
    #
    #         val = p.D*(1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
    #         col = row
    #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
    #     end
    # end
    # use Robin BC
    for vi in 1:p.num_freq
        for j in 1:p.n_rot
            row = p.layer_unknown*Ri + (vi-1)*p.n_rot + j
            s = put_row_col_val(rowind, colind, value, row, row, -p.D/p.Δr, s)
            s = put_row_col_val(rowind, colind, value, row, row-p.layer_unknown, p.D/p.Δr, s)
        end
    end

end
