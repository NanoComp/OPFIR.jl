function put_row_col_val(rowind, colind, value, row, col, val, s)
    rowind[s] = row
    colind[s] = col
    value[s] = val
    return (s+1)
end

function compute_row_col_val(rowind, colind, value, p, sol_0)
    s = 1
    # rotational levels in V0 and V3
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
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

            # col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1 #V0A
            # val = - p.pump_IR[vi, ri] * p.CL * p.gauss_dist[vi]
            # s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
            val = p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            # col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2 #V3A
            # val = p.pump_IR[vi, ri] * p.CU * p.gauss_dist[vi] * p.g_L/p.g_U
            # s = put_row_col_val(rowind, colind, value, row, col, val, s)

            ### pumping for U level
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+1) #
            val = p.pump_IR[vi, ri]
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
            #
            # col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1 # V0A
            # val = p.pump_IR[vi, ri] * p.CL * p.gauss_dist[vi]
            # s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = row
            val = - p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
            #
            # col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2 # V3A
            # val = - p.pump_IR[vi, ri] * p.CU * p.gauss_dist[vi] * p.g_L/p.g_U
            # s = put_row_col_val(rowind, colind, value, row, col, val, s)
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
        # rotational levels
        index_diffu = vcat(p.layer_unknown * (ri-1) + 1 : p.layer_unknown * ri - 1)
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
    index_diffu = vcat(1:p.layer_unknown-1)
    for k in index_diffu
        row = k
        val = p.D*(-1.0/(p.Δr)^2 - 1.0/2.0/p.Δr/p.r_int[1])
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        col = row + p.layer_unknown
        val = p.D*(1.0/(p.Δr)^2 + 1.0/2/p.Δr/p.r_int[1])
        s = put_row_col_val(rowind, colind, value, row, col, val, s)
    end

    ################ question: BC for the wall? #######
    Ri = p.num_layers
    # Robin BC:
    vbar = p.v_avg / sqrt(2) / p.norm_time/2
    # rotatinal levels
    for row in vcat(p.layer_unknown*Ri+1:p.layer_unknown*(Ri+1)-1] # index offset, or starting index
        val = -p.D/p.Δr + vbar * (1-p.f)/4
        s = put_row_col_val(rowind, colind, value, row, row-p.layer_unknown, val, s)

        val = p.D/p.Δr + vbar * (1-p.f_G_0)/4
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = -vbar*p.f_G_0/4
        for k in [1, 2]
            s = put_row_col_val(rowind, colind, value, row, row-p.n_vib+k, val, s)
            s = put_row_col_val(rowind, colind, value, row, row+k, val, s)
        end
    end

    # V3A and V3E at N+1 grid
    for row in [p.layer_unknown * Ri + 2, p.layer_unknown * Ri + p.n_vib÷2 + 2] # index offset, or starting index
        val = -p.D/p.Δr + vbar * (1-p.f_3_0)/4
        s = put_row_col_val(rowind, colind, value, row, row-p.n_vib, val, s)

        val = p.D/p.Δr + vbar * (1-p.f_3_0)/4
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = -vbar*p.f_3_0/4
        for k in [-1, 1]
            s = put_row_col_val(rowind, colind, value, row, row-p.n_vib+k, val, s)
            s = put_row_col_val(rowind, colind, value, row, row+k, val, s)
        end
    end
    # VΣA and VΣE at N+1 grid
    for row in [p.layer_unknown * Ri + 3, p.layer_unknown * Ri + p.n_vib÷2 + 3] # index offset, or starting index
        val = -p.D/p.Δr + vbar * (1-p.f_6_0)/4
        s = put_row_col_val(rowind, colind, value, row, row-p.n_vib, val, s)

        val = p.D/p.Δr + vbar * (1-p.f_6_0)/4
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = -vbar*p.f_6_0/4
        for k in [-1, -2]
            s = put_row_col_val(rowind, colind, value, row, row-p.n_vib+k, val, s)
            s = put_row_col_val(rowind, colind, value, row, row+k, val, s)
        end
    end


    if p.model_flag==2
        ### for case gamma -> infinity
        for ri in 1:p.num_layers
            # V3A and V3E:
            rowA = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
            for row in vcat(rowA, rowA + p.n_vib÷2)
                s = put_row_col_val(rowind, colind, value, row, row-1, 1., s)
                s = put_row_col_val(rowind, colind, value, row, row,   1., s)
                s = put_row_col_val(rowind, colind, value, row, row+1, 1., s)
            end
            # VΣA and VΣE
            rowA = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 3
            rowE = rowA + p.n_vib÷2
            for row in vcat(rowA, rowE)
                s = put_row_col_val(rowind, colind, value, row, row, 1., s)
            end
            s = put_row_col_val(rowind, colind, value, rowA, rowA-1, -p.k36A[ri]/p.k63A[ri], s)
            s = put_row_col_val(rowind, colind, value, rowE, rowE-1, -p.k36E[ri]/p.k63E[ri], s)
        end
    end

end
