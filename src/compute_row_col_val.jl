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
            for j in 1:p.n_rot
                row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + j
                val = -p.ka[j] - p.kro
                s = put_row_col_val(rowind, colind, value, row, row, val, s)
            end

            for li in 1:p.n_rot
                for lj in 1:p.n_rot
                    if p.kDDmat[li, lj] > 0
                        row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + li
                        col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + lj
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
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-2)
            val = - p.pump_IR[vi, ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1 #V0
            val = - p.pump_IR[vi, ri] * p.CL * p.gauss_dist[vi]
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-2)
            val = p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2 #V3
            val = p.pump_IR[vi, ri] * p.CU * p.gauss_dist[vi] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            ### pumping for U level
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-2)
            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-2) #
            val = p.pump_IR[vi, ri]
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1 # V0
            val = p.pump_IR[vi, ri] * p.CL * p.gauss_dist[vi]
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = row
            val = - p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2 # V3
            val = - p.pump_IR[vi, ri] * p.CU * p.gauss_dist[vi] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
        end
    end

    ### stimulated emission term:
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-2) #V0, J4 level
            val = - p.WiL*(p.g_L+2)/p.g_L ## J4 -> J5 (diagonal):
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = p.WiL ## J5 -> J4:
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-1) #V0, J5 level
            val = - p.WiL  # J5 -> J4 (diagonal)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = p.WiL*(p.g_L+2)/p.g_L # J4 -> J5
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-3) #V3, J4
            val = - p.WiU*p.g_U/(p.g_U-2) ## J4 -> J5 (diagonal)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = p.WiU ## J5 -> J4:
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-2) #V3, J5
            val = - p.WiU ## J5 -> J4
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = p.WiU*p.g_U/(p.g_U-2) ## J4 -> J5
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end
    end

    #### Thermal pools: ##########
    ### from rotational levels to thermal pools ###
    for ri in 1:p.num_layers
        # V0 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
        for vi in 1:p.num_freq
            for k in 1:p.n_rot÷2
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + k
                val = p.ka[k] + p.kro/2
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
        end

        # V0 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + p.n_vib÷2 + 1
        for vi in 1:p.num_freq
            for k in 1:p.n_rot÷2
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + k
                val = p.kro/2
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
        end

        if p.model_flag == 1
            # V3 A type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
            for vi in 1:p.num_freq
                for k in 1:p.n_rot÷2
                    col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + k
                    val = p.ka[k+p.n_rot÷2] + p.kro/2
                    s = put_row_col_val(rowind, colind, value, row, col, val, s)
                end
            end

            # V3 E type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + p.n_vib÷2 + 2
            for vi in 1:p.num_freq
                for k in 1:p.n_rot÷2
                    col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + k
                    val = p.kro/2
                    s = put_row_col_val(rowind, colind, value, row, col, val, s)
                end
            end
        end
    end

    ### transitions between thermal pools ###
    ## A type:
    for ri in 1:p.num_layers
        rot_idx = 1
        # V0 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1

        # V3 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1
        #### transtion between V3 and V6 ####
        if p.model_flag==1 # || p.model_flag==2
            val = -p.k36A[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k63A[ri]
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        elseif p.model_flag==2 #&& p.solstart_flag==0
            # val = -p.k36A[ri]
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            # val = +p.k63A[ri]
            # s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        end

        # V6 or V_Σ A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1
        #### transtion between V3 and V6: ####
        if p.model_flag==1
            val = -p.k63A[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k36A[ri]
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        elseif p.model_flag==2 #&& p.solstart_flag==0
            # val = -p.k63A[ri]
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            # val = +p.k36A[ri]
            # s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end

        if p.model_flag==1
            # V23 A type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
            rot_idx += 1
            #### transitions between V23 and V36 ####
            val = -p.k2336
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k3623
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # V36 A type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
            rot_idx += 1
            #### transition between V36 and V23, V26 ####
            val = -p.k3623 - p.k3626
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k2336
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
            val = +p.k2636
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # V26 A type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
            rot_idx += 1
            #### transition between V26 and V36: ####
            val = -p.k2636
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k3626
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end

################################################################################
        # V0 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1

        # V3 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1
        #### transtion between V3 and V6 ####
        if p.model_flag==1
            val = -p.k36E[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k63E[ri]
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        elseif p.model_flag==2
            # val = -p.k36E[ri]
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            # val = +p.k63E[ri]
            # s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        end

        # V6 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1
        #### transtion between V3 and V6: ####
        if p.model_flag==1
            val = -p.k63E[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k36E[ri]
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        elseif p.model_flag==2
            # val = -p.k63E[ri]
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            # val = +p.k36E[ri]
            # s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end

        if p.model_flag==1
            # V23 E type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
            rot_idx += 1
            #### transitions between V23 and V36 ####
            val = -p.k2336
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k3623
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # V36 E type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
            rot_idx += 1
            #### transition between V36 and V23, V26 ####
            val = -p.k3623 - p.k3626
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k2336
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
            val = +p.k2636
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # V26 E type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
            rot_idx += 1
            #### transition between V26 and V36: ####
            val = -p.k2636
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k3626
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end # if model_flag==1
    end # for ri, thermal pools

    ########################## add the diffusion term ##########################

    for ri in 2:p.num_layers-1
        if p.model_flag==2
            index_diffu = vcat(p.layer_unknown * (ri-1) + 1 : p.layer_unknown * ri - 5,
                               p.layer_unknown*ri-2)
        else
            index_diffu = vcat(p.layer_unknown * (ri-1) + 1 : p.layer_unknown * ri - 7,
                               p.layer_unknown*ri-5:p.layer_unknown*ri-1)
        end
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
    if p.model_flag==2
        index_diffu = vcat(1:p.layer_unknown-5, p.layer_unknown-2)
    else
        index_diffu = vcat(1:p.layer_unknown-7,p.layer_unknown-5:p.layer_unknown-1)
    end
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
    # for rotational levels, just use Neumann BC
    index_diffu = p.layer_unknown * (Ri-1) + 1 : p.layer_unknown * Ri - p.n_vib
    for k in index_diffu
        row = k
        col = row - p.layer_unknown
        val = p.D*(-1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
        s = put_row_col_val(rowind, colind, value, row, col, val, s)

        val = p.D*(-2.0/(p.Δr)^2) + p.D*(1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
    end
    # for vibrational levels, use Robin BC
    if p.model_flag == 2
        for k in vcat(p.layer_unknown * Ri - p.n_vib + 1, p.layer_unknown * Ri - 2)
            row = k
            col = row - p.layer_unknown
            val = p.D*(-1.0/(2*p.Δr*p.r_int[Ri]) + 1.0/(p.Δr)^2)
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            val = p.D*(-2.0/(p.Δr)^2)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = p.D*(1.0/(2*p.Δr*p.r_int[Ri]) + 1.0/(p.Δr)^2)
            col = row + p.n_vib #!!!
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
        end
        # Robin BC:
        vbar = p.v_avg / sqrt(2) / p.norm_time/2
        # V0A and V0E at N+1 grid
        for row in [p.layer_unknown * Ri + 1, p.layer_unknown * Ri + p.n_vib÷2 + 1] # index offset, or starting index
            val = -p.D/p.Δr + vbar * (1-p.f_G_0)/4
            s = put_row_col_val(rowind, colind, value, row, row-p.n_vib, val, s)

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
    else
        index_diffu = vcat(p.layer_unknown * Ri - p.n_vib + 1:p.layer_unknown * Ri - p.n_vib + 5,
        p.layer_unknown * Ri - 5: p.layer_unknown * Ri - 1)
        for k in index_diffu
            row = k
            col = row - p.layer_unknown
            val = p.D*(-1.0/(2*p.Δr*p.r_int[Ri]) + 1.0/(p.Δr)^2)
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            val = p.D*(-2.0/(p.Δr)^2)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = p.D*(1.0/(2*p.Δr*p.r_int[Ri]) + 1.0/(p.Δr)^2)
            col = row + p.n_vib #!!!
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
        end
        # Robin BC:
        vbar = p.v_avg / sqrt(2) / p.norm_time/2
        # V0A and V0E at N+1 grid
        for row in [p.layer_unknown * Ri + 1, p.layer_unknown * Ri + p.n_vib÷2 + 1] # index offset, or starting index
            val = -p.D/p.Δr + vbar * (1-p.f_G_0)/4
            s = put_row_col_val(rowind, colind, value, row, row-p.n_vib, val, s)

            val = p.D/p.Δr + vbar * (1-p.f_G_0)/4
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = -vbar*p.f_G_0/4
            for k in [1, 2, 3, 4, 5]
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
            for k in [-1, 1, 2, 3, 4]
                s = put_row_col_val(rowind, colind, value, row, row-p.n_vib+k, val, s)
                s = put_row_col_val(rowind, colind, value, row, row+k, val, s)
            end
        end
        # V6A and V6E at N+1 grid
        for row in [p.layer_unknown * Ri + 3, p.layer_unknown * Ri + p.n_vib÷2 + 3] # index offset, or starting index
            val = -p.D/p.Δr + vbar * (1-p.f_6_0)/4
            s = put_row_col_val(rowind, colind, value, row, row-p.n_vib, val, s)

            val = p.D/p.Δr + vbar * (1-p.f_6_0)/4
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = -vbar*p.f_6_0/4
            for k in [-2, -1, 1, 2, 3]
                s = put_row_col_val(rowind, colind, value, row, row-p.n_vib+k, val, s)
                s = put_row_col_val(rowind, colind, value, row, row+k, val, s)
            end
        end
        # 2V3A and 2V3E at N+1 grid
        for row in [p.layer_unknown * Ri + 4, p.layer_unknown * Ri + p.n_vib÷2 + 4] # index offset, or starting index
            val = -p.D/p.Δr + vbar * (1-p.f_23)/4
            s = put_row_col_val(rowind, colind, value, row, row-p.n_vib, val, s)

            val = p.D/p.Δr + vbar * (1-p.f_23)/4
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = -vbar*p.f_23/4
            for k in [-3, -2, -1, 1, 2]
                s = put_row_col_val(rowind, colind, value, row, row-p.n_vib+k, val, s)
                s = put_row_col_val(rowind, colind, value, row, row+k, val, s)
            end
        end
        # V36A and V36E at N+1 grid
        for row in [p.layer_unknown * Ri + 5, p.layer_unknown * Ri + p.n_vib÷2 + 5] # index offset, or starting index
            val = -p.D/p.Δr + vbar * (1-p.f_36)/4
            s = put_row_col_val(rowind, colind, value, row, row-p.n_vib, val, s)

            val = p.D/p.Δr + vbar * (1-p.f_36)/4
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = -vbar*p.f_36/4
            for k in [-4, -3, -2, -1, 1]
                s = put_row_col_val(rowind, colind, value, row, row-p.n_vib+k, val, s)
                s = put_row_col_val(rowind, colind, value, row, row+k, val, s)
            end
        end
        # V26A and V26E at N+1 grid
        for row in [p.layer_unknown * Ri + 6, p.layer_unknown * Ri + p.n_vib÷2 + 6] # index offset, or starting index
            # val = -p.D/p.Δr + vbar * (1-p.f_26)/4
            # s = put_row_col_val(rowind, colind, value, row, row-p.n_vib, val, s)
            #
            # val = p.D/p.Δr + vbar * (1-p.f_26)/4
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            #
            # val = -vbar*p.f_26/4
            for k in [-5, -4, -3, -2, -1, 0]
                # s = put_row_col_val(rowind, colind, value, row, row-p.n_vib+k, val, s)
                s = put_row_col_val(rowind, colind, value, row, row+k, 1., s)
            end
        end
    end

    ############## add the V-swap process V0A + V3E <-> V_3A + V0E> ##########
    N0A_0 = p.ntotal * p.f_G_0/2
    N0E_0 = N0A_0
    N3A_0 = p.ntotal * p.f_3_0/2
    N3E_0 = N3A_0
    non_kvs = p.kvs
    # non_kvs = 0
    for ri in 1:p.num_layers
        offset = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot
        N0A = sol_0[offset + 1]
        N3A = sol_0[offset + 2]
        N0E = sol_0[offset + p.n_vib÷2 + 1]
        N3E = sol_0[offset + p.n_vib÷2 + 2]

        ## V0A:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1

        val = -non_kvs * N3E - non_kvs * N3E_0
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = -non_kvs * N0A_0
        s = put_row_col_val(rowind, colind, value, row, row+p.n_vib÷2+1, val, s)

        val = + non_kvs * N3A + non_kvs * N3A_0
        s = put_row_col_val(rowind, colind, value, row, row+p.n_vib÷2, val, s)

        val = + non_kvs * N0E_0
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        if p.model_flag == 1 # for model 2, V3A and V3E has conservation eq
            # ## V3E
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + p.n_vib÷2 + 2

            val = -non_kvs * N0A - non_kvs * N0A_0
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = -non_kvs * N3E_0
            s = put_row_col_val(rowind, colind, value, row, row-p.n_vib÷2-1, val, s)

            val = + non_kvs * N3A_0
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            val = + non_kvs * N0E + non_kvs * N0E_0
            s = put_row_col_val(rowind, colind, value, row, row-p.n_vib÷2, val, s)

            ## V3A:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2

            val = -non_kvs * N0E - non_kvs * N0E_0
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = -non_kvs * N3A_0
            s = put_row_col_val(rowind, colind, value, row, row+p.n_vib÷2-1, val, s)

            val = + non_kvs * N0A + non_kvs * N0A_0
            s = put_row_col_val(rowind, colind, value, row, row+p.n_vib÷2, val, s)

            val = + non_kvs * N3E_0
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end

        ## V0E:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + p.n_vib÷2 + 1
        val = -non_kvs * N3A - non_kvs * N3A_0
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = -non_kvs * N0E_0
        s = put_row_col_val(rowind, colind, value, row, row-p.n_vib÷2+1, val, s)

        val = +non_kvs * N0A_0
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        val = + non_kvs * N3E + non_kvs * N3E_0
        s = put_row_col_val(rowind, colind, value, row, row-p.n_vib÷2, val, s)

        if p.model_flag == 1
            ## V-split process: V0 + V23 <-> V3 + V3 for 6-level model
            non_kvsplit = p.kvs / p.σ_VS * 12.4 * 0
            N23_0 = p.ntotal * p.f_23/2
            N23A = sol_0[offset + 4]
            N23E = sol_0[offset + p.n_vib÷2+ 4]

            # V0A:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1

            val = -non_kvsplit * (N23A + N23_0)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = -non_kvsplit * N0A_0
            s = put_row_col_val(rowind, colind, value, row, row+3, val, s)

            val = + non_kvsplit * (N3A + 2 * N3A_0)
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # V23A:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 4

            val = -non_kvsplit * (N0A + N0A_0)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = -non_kvsplit * N23_0
            s = put_row_col_val(rowind, colind, value, row, row-3, val, s)

            val = + non_kvsplit * (N3A + 2 * N3A_0)
            s = put_row_col_val(rowind, colind, value, row, row-2, val, s)

            # V3A:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2

            val = -non_kvsplit * (N23A + N23_0)
            s = put_row_col_val(rowind, colind, value, row, row-1, -2.*val, s)

            val = -non_kvsplit * N0A_0
            s = put_row_col_val(rowind, colind, value, row, row+2, -2.*val, s)

            val = + non_kvsplit * (N3A + 2 * N3A_0)
            s = put_row_col_val(rowind, colind, value, row, row, -2.*val, s)

            # V0E:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + p.n_vib÷2 + 1

            val = -non_kvsplit * (N23E + N23_0)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = -non_kvsplit * N0E_0
            s = put_row_col_val(rowind, colind, value, row, row+3, val, s)

            val = + non_kvsplit * (N3E + 2 * N3E_0)
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # V23E:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + p.n_vib÷2 + 4

            val = -non_kvsplit * (N0E + N0E_0)
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = -non_kvsplit * N23_0
            s = put_row_col_val(rowind, colind, value, row, row-3, val, s)

            val = + non_kvsplit * (N3E + 2 * N3E_0)
            s = put_row_col_val(rowind, colind, value, row, row-2, val, s)

            # V3E:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + p.n_vib÷2 + 2

            val = -non_kvsplit * (N23E + N23_0)
            s = put_row_col_val(rowind, colind, value, row, row-1, -2.*val, s)

            val = -non_kvsplit * N0E_0
            s = put_row_col_val(rowind, colind, value, row, row+2, -2.*val, s)

            val = + non_kvsplit * (N3E + 2 * N3E_0)
            s = put_row_col_val(rowind, colind, value, row, row, -2.*val, s)
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
