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
            # oscil: ground vib: J3
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 1
            ## J3 -> J4, and thermal pool (SPT and VS) (diagonal):
            val = -p.k12_G - p.k1a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J4 -> J3:
            val = p.k21_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # oscil: ground vib: J4
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 2
            ## J4 -> J5, J3 (diagonal):
            val = - p.k23_G - p.k21_G - p.k2a - p.kro - p.WiL*p.g_U/p.g_L
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J5, J3 -> J4:
            val = p.k32_G + p.WiL
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k12_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
            ### pumping:
            val = -p.pump_IR[vi, ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = - p.pump_IR[vi, ri] * p.C4L * p.gauss_dist[vi]
            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + 3
            val = p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
            val = p.pump_IR[vi, ri] * p.C5U * p.gauss_dist[vi] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            # oscil: ground vib: J5
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 3
            ## J5 -> J6, J4 (diagonal):
            val = - p.k34_G - p.k32_G - p.k3a - p.kro - p.WiL
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J6, J4 -> J5:
            val = p.k43_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k23_G + p.WiL*p.g_U/p.g_L
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J6
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 4
            ## J6 -> J7, J5 (diagonal):
            val = - p.k45_G - p.k43_G - p.k4a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J7, J5 -> J6:
            val = p.k54_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k34_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J7
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 5
            ## J7 -> J8, J6 (diagonal):
            val = - p.k56_G - p.k54_G - p.k5a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J8, J6 -> J7:
            val = p.k65_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k45_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J8
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 6
            ## J8 -> J9, J7 (diagonal):
            val = - p.k67_G - p.k65_G - p.k6a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J9, J7 -> J8:
            val = p.k76_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k56_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J9
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 7
            ## J9 -> J10, J8 (diagonal):
            val = - p.k78_G - p.k76_G - p.k7a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J10, J8 -> J9:
            val = p.k87_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k67_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J10
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 8
            ## J10 -> J11, J9 (diagonal):
            val = - p.k89_G - p.k87_G - p.k8a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J11, J9 -> J11:
            val = p.k98_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k78_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J11
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9
            ## J11 -> J10 (diagonal):
            val = - p.k98_G - p.k9a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J10 -> J11:
            val = p.k89_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end
    end
    # rotational levels in V3
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
            # oscil: 3 vib: J3
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 1
            ## J3 -> J4 (diagonal):
            val = -p.k12_3 - p.k10a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J4 -> J3:
            val = p.k21_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # oscil: 3 vib: J4
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 2
            ## J4 -> J5, J3 (diagonal):
            val = - p.k23_3 - p.k21_3 - p.k11a - p.kro - p.WiU*p.g_U/p.g_L
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J5, J3 -> J4:
            val = p.k32_3 + p.WiU
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k12_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J5
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 3
            ## J5 -> J6, J4 (diagonal):
            val = - p.k34_3 - p.k32_3 - p.k12a - p.kro - p.WiU
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J6, J4 -> J5:
            val = p.k43_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k23_3 + p.WiU*p.g_U/p.g_L
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            ### pumping:
            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 2
            val = p.pump_IR[vi, ri]
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
            val = p.pump_IR[vi, ri] * p.C4L * p.gauss_dist[vi]
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = row
            val = - p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
            val = - p.pump_IR[vi, ri] * p.C5L * p.gauss_dist[vi] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            # vib: 3 rot: J6
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 4
            ## J6 -> J7, J5 (diagonal):
            val = - p.k45_3 - p.k43_3 - p.k13a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J7, J5 -> J6:
            val = p.k54_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k34_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J7
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 5
            ## J7 -> J8, J6 (diagonal):
            val = - p.k56_3 - p.k54_3 - p.k14a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J8, J6 -> J7:
            val = p.k65_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k45_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J8
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 6
            ## J8 -> J9, J7 (diagonal):
            val = - p.k67_3 - p.k65_3 - p.k15a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J9, J7 -> J8:
            val = p.k76_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k56_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J9
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 7
            ## J9 -> J10, J8 (diagonal):
            val = - p.k78_3 - p.k76_3 - p.k16a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J10, J8 -> J9:
            val = p.k87_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k67_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J10
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 8
            ## J10 -> J11, J9 (diagonal):
            val = - p.k89_3 - p.k87_3 - p.k17a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J11, J9 -> J11:
            val = p.k98_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k78_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J11
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 9
            ## J11 -> J10 (diagonal):
            val = - p.k98_3 - p.k18a - p.kro
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J10 -> J11:
            val = p.k89_3
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
                val = p.k1a + p.kro/2
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
        end

        # V3 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
        for vi in 1:p.num_freq
            for k in 1:p.n_rot÷2
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + k
                val = p.k1a + p.kro/2
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

    ### transitions between thermal pools ###
    ## A type:
    for ri in 1:p.num_layers
        rot_idx = 1
        # V0 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1
        #
        # val =  -p.kwall[ri] + p.f_G_0 * p.kwall[ri]
        # s = put_row_col_val(rowind, colind, value, row, row, val, s)
        #
        # val = p.f_G_0 * p.kwall[ri]
        # if p.model_flag==1
        #     for col in [row+1; row+2; row+3; row+4; row+5]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # elseif p.model_flag==2
        #     for col in [row+1; row+2]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # end
        #
        # V3 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1
        #
        # val = - p.kwall[ri] + p.f_3_0 * p.kwall[ri]
        # s = put_row_col_val(rowind, colind, value, row, row, val, s)
        #
        # val = p.f_3_0 * p.kwall[ri]
        # if p.model_flag==1
        #     for col in [row-1; row+1; row+2; row+3; row+4]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # elseif p.model_flag==2
        #     for col in [row-1; row+1]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # end

        #### transtion between V3 and V6 ####
        if p.model_flag==1 # || p.model_flag==2
            val = -p.k36A[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k63A[ri]
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        elseif p.model_flag==2 #&& p.solstart_flag==0
            val = -p.k36A[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k63A[ri]
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        end
        # V6 or V_Σ A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1
        #
        # val = -p.kwall[ri] + p.f_6_0 * p.kwall[ri]
        # s = put_row_col_val(rowind, colind, value, row, row, val, s)
        #
        # val = p.f_6_0 * p.kwall[ri]
        # if p.model_flag==1
        #     for col in [row-2; row-1; row+1; row+2; row+3]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # elseif p.model_flag==2
        #     for col in [row-2; row-1]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # end

        #### transtion between V3 and V6: ####
        if p.model_flag==1
            val = -p.k63A[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k36A[ri]
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        elseif p.model_flag==2 #&& p.solstart_flag==0
            val = -p.k63A[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k36A[ri]
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end

        if p.model_flag==1
            # V23 A type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
            rot_idx += 1
            #
            # val = -p.kwall[ri] + p.f_23 * p.kwall[ri]
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            #
            # val = p.f_23 * p.kwall[ri]
            # for col in [row-3; row-2; row-1; row+1; row+2]
            #     s = put_row_col_val(rowind, colind, value, row, col, val, s)
            # end
            #### transitions between V23 and V36 ####
            val = -p.k2336
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k3623
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # V36 A type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
            rot_idx += 1

            # val = -p.kwall[ri] + p.f_36 * p.kwall[ri]
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            #
            # val = p.f_36 * p.kwall[ri]
            # for col in [row-4; row-3; row-2; row-1; row+1]
            #     s = put_row_col_val(rowind, colind, value, row, col, val, s)
            # end

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

            # val = -p.kwall[ri] + p.f_26 * p.kwall[ri]
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            #
            # val = p.f_26 * p.kwall[ri]
            # for col in [row-5; row-4; row-3; row-2; row-1]
            #     s = put_row_col_val(rowind, colind, value, row, col, val, s)
            # end

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
        #
        # val = -p.kwall[ri] + p.f_G_0 * p.kwall[ri]
        # s = put_row_col_val(rowind, colind, value, row, row, val, s)
        #
        # val = p.f_G_0 * p.kwall[ri]
        # if p.model_flag==1
        #     for col in [row+1; row+2; row+3; row+4; row+5]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # elseif p.model_flag==2
        #     for col in [row+1; row+2]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # end

        # V3 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1
        #
        # val = - p.kwall[ri] + p.f_3_0 * p.kwall[ri]
        # s = put_row_col_val(rowind, colind, value, row, row, val, s)
        #
        # val = p.f_3_0 * p.kwall[ri]
        # if p.model_flag==1
        #     for col in [row-1; row+1; row+2; row+3; row+4]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # elseif p.model_flag==2
        #     for col in [row-1; row+1]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # end

        #### transtion between V3 and V6 ####
        if p.model_flag==1
            val = -p.k36E[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k63E[ri]
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        elseif p.model_flag==2
            # put in the rhs
            val = -p.k36E[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k63E[ri]
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        end

        # V6 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
        rot_idx += 1
        #
        # val = -p.kwall[ri] + p.f_6_0 * p.kwall[ri]
        # s = put_row_col_val(rowind, colind, value, row, row, val, s)
        #
        # val = p.f_6_0 * p.kwall[ri]
        # if p.model_flag==1
        #     for col in [row-2; row-1; row+1; row+2; row+3]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # elseif p.model_flag==2
        #     for col in [row-2; row-1]
        #         s = put_row_col_val(rowind, colind, value, row, col, val, s)
        #     end
        # end

        #### transtion between V3 and V6: ####
        if p.model_flag==1
            val = -p.k63E[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k36E[ri]
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        elseif p.model_flag==2
            # put it in the rhs
            val = -p.k63E[ri]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k36E[ri]
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end

        if p.model_flag==1
            # V23 E type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
            rot_idx += 1
            #
            # val = -p.kwall[ri] + p.f_23 * p.kwall[ri]
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            #
            # val = p.f_23 * p.kwall[ri]
            # for col in [row-3; row-2; row-1; row+1; row+2]
            #     s = put_row_col_val(rowind, colind, value, row, col, val, s)
            # end
            #### transitions between V23 and V36 ####
            val = -p.k2336
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k3623
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # V36 E type:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + rot_idx
            rot_idx += 1
            #
            # val = -p.kwall[ri] + p.f_36 * p.kwall[ri]
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            #
            # val = p.f_36 * p.kwall[ri]
            # for col in [row-4; row-3; row-2; row-1; row+1]
            #     s = put_row_col_val(rowind, colind, value, row, col, val, s)
            # end

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
            #
            # val = -p.kwall[ri] + p.f_26 * p.kwall[ri]
            # s = put_row_col_val(rowind, colind, value, row, row, val, s)
            #
            # val = p.f_26 * p.kwall[ri]
            # for col in [row-5; row-4; row-3; row-2; row-1]
            #     s = put_row_col_val(rowind, colind, value, row, col, val, s)
            # end

            #### transition between V26 and V36: ####
            val = -p.k2636
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            val = +p.k3626
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end

        # println(rot_idx)
    end

    ########################## add the diffusion term ##########################

    for ri in 2:p.num_layers-1
        index_diffu = p.layer_unknown * (ri-1) + 1 : p.layer_unknown * ri
        # dx1 = p.r_int[ri] - p.r_int[ri-1]
        # dx2 = p.r_int[ri+1] - p.r_int[ri]
        # ai1 = -dx2/dx1/(dx1+dx2)
        # bi1 = (dx2-dx1)/dx1/dx2
        # ci1 = dx1/dx2/(dx1+dx2)
        # ai2 = 2/dx1/(dx1+dx2)
        # bi2 = -2/dx1/dx2
        # ci2 = 2/dx2/(dx1+dx2)
        for k in index_diffu
            row = k
            col = row - p.layer_unknown
            val = p.D*(-1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
            # val = p.D * (ai2 + ai1/p.r_int[ri])
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            val = p.D*(-2.0/(p.Δr)^2)
            # val = p.D * (bi2 + bi1/p.r_int[ri])
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = p.D*(1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
            # val = p.D * (ci2 + ci1/p.r_int[ri])
            col = row + p.layer_unknown
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
        end
    end

    ## first layer
    # dx1 = p.r_int[1] * 2
    # dx2 = p.r_int[2] - p.r_int[1]
    # ai1 = -dx2/dx1/(dx1+dx2)
    # bi1 = (dx2-dx1)/dx1/dx2
    # ci1 = dx1/dx2/(dx1+dx2)
    # ai2 = 2/dx1/(dx1+dx2)
    # bi2 = -2/dx1/dx2
    # ci2 = 2/dx2/(dx1+dx2)
    for k in 1:p.layer_unknown
        row = k
        val = p.D*(-1.0/(p.Δr)^2 - 1.0/2.0/p.Δr/p.r_int[1])
        # val = p.D * (ai2+bi2 + ai1/p.r_int[1]+bi1/p.r_int[1])
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        col = row + p.layer_unknown
        val = p.D*(1.0/(p.Δr)^2 + 1.0/2/p.Δr/p.r_int[1])
        # val = p.D * (ci2 + ci1/p.r_int[1])
        s = put_row_col_val(rowind, colind, value, row, col, val, s)
    end

    ################ question: BC for the wall? #######
    ri = p.num_layers
    # for rotational levels, just use Neumann BC
    index_diffu = p.layer_unknown * (ri-1) + 1 : p.layer_unknown * ri - p.n_vib
    for k in index_diffu
        row = k
        col = row - p.layer_unknown
        val = p.D*(-1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
        # val = p.D * (ai2 + ai1/p.r_int[end])
        s = put_row_col_val(rowind, colind, value, row, col, val, s)

        val = p.D*(-2.0/(p.Δr)^2) + p.D*(1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
        # val = p.D * (bi2+ci2 + bi1/p.r_int[end] + ci1/p.r_int[end])
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
    end
    # for vibrational levels, use Robin BC
    for k in p.layer_unknown * ri - p.n_vib + 1 : p.layer_unknown * ri
        row = k
        col = row - p.layer_unknown
        val = p.D*(-1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
        s = put_row_col_val(rowind, colind, value, row, col, val, s)

        val = p.D*(-2.0/(p.Δr)^2)
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.D*(1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
        col = row + p.n_vib #!!!
        s = put_row_col_val(rowind, colind, value, row, col, val, s)
    end
    # Robin BC:
    vbar = p.v_avg / sqrt(2) / p.norm_time/2
    # V0A and V0E at N+1 grid
    for row in [p.layer_unknown * ri + 1, p.layer_unknown * ri + p.n_vib÷2 + 1] # index offset, or starting index
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
    for row in [p.layer_unknown * ri + 2, p.layer_unknown * ri + p.n_vib÷2 + 2] # index offset, or starting index
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
    for row in [p.layer_unknown * ri + 3, p.layer_unknown * ri + p.n_vib÷2 + 3] # index offset, or starting index
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

    #
    #
    # dx1 = p.r_int[end] - p.r_int[end-1]
    # dx2 = (p.radius/100 - p.r_int[end]) * 2
    # ai1 = -dx2/dx1/(dx1+dx2)
    # bi1 = (dx2-dx1)/dx1/dx2
    # ci1 = dx1/dx2/(dx1+dx2)
    # ai2 = 2/dx1/(dx1+dx2)
    # bi2 = -2/dx1/dx2
    # ci2 = 2/dx2/(dx1+dx2)
    # for k in p.layer_unknown*(p.num_layers-1)+1:p.layer_unknown*p.num_layers
    #     row = k
    #     col = row - p.layer_unknown
    #     # val = p.D*(-1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
    #     val = p.D * (ai2 + ai1/p.r_int[end])
    #     s = put_row_col_val(rowind, colind, value, row, col, val, s)
    #
    #     # val = p.D*(-2.0/(p.Δr)^2) + p.D*(1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
    #     val = p.D * (bi2+ci2 + bi1/p.r_int[end] + ci1/p.r_int[end])
    #     s = put_row_col_val(rowind, colind, value, row, row, val, s)
    # end

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
        ## non-thermal rotational levels:
        for vi in 1:p.num_freq
            for ell in 1:p.n_rot÷2 # in V0A
                row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + ell
                val = -non_kvs * (N3E_0 + N3E)
                s = put_row_col_val(rowind, colind, value, row, row, val, s)
            end
            for ell in 1:p.n_rot÷2 # in V3A
                row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + ell
                val = -non_kvs * (N0E_0 + N0E)
                s = put_row_col_val(rowind, colind, value, row, row, val, s)
            end
        end

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

        # from non-thermal:
        for vi in 1:p.num_freq
            for ell in 1:p.n_rot÷2 # in V3A
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + ell
                val = non_kvs * (N0E_0 + N0E)
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
        end

        ## V3E
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + p.n_vib÷2 + 2

        val = -non_kvs * N0A - non_kvs * N0A_0
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = -non_kvs * N3E_0
        s = put_row_col_val(rowind, colind, value, row, row-p.n_vib÷2-1, val, s)

        val = + non_kvs * N3A_0
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

        val = + non_kvs * N0E + non_kvs * N0E_0
        s = put_row_col_val(rowind, colind, value, row, row-p.n_vib÷2, val, s)

        # from non-thermal
        for vi in 1:p.num_freq
            for ell in 1:p.n_rot÷2 # in V0A
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + ell
                val = -non_kvs * (N3E_0 + N3E)
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
            for ell in 1:p.n_rot÷2 # in V3A
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + ell
                val = non_kvs * (N0E_0 + N0E)
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
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

        # from non-thermal
        for vi in 1:p.num_freq
            for ell in 1:p.n_rot÷2 # in V0A
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + ell
                val = non_kvs * (N3E_0 + N3E)
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
            for ell in 1:p.n_rot÷2 # in V3A
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + ell
                val = -non_kvs * (N0E_0 + N0E)
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
        end

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

        # from non-thermal
        for vi in 1:p.num_freq
            for ell in 1:p.n_rot÷2 # in V0A
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + ell
                val = non_kvs * (N3E_0 + N3E)
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
        end


        if p.model_flag == 1
            ## V-split process: V0 + V23 <-> V3 + V3 for 6-level model
            non_kvsplit = p.kvs / p.σ_VS * 12.4
            N23_0 = p.ntotal * p.f_23/2

            # N0A = sol_0[offset + 1]
            N23A = sol_0[offset + 4]
            # N3A = sol_0[offset + 2]
            # N0E = sol_0[offset + p.n_vib÷2+ 1]
            N23E = sol_0[offset + p.n_vib÷2+ 4]
            # N3E = sol_0[offset + p.n_vib÷2+ 2]
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
            s = put_row_col_val(rowind, colind, -2.*value, row, row-1, val, s)

            val = -non_kvsplit * N0A_0
            s = put_row_col_val(rowind, colind, -2.*value, row, row+2, val, s)

            val = + non_kvsplit * (N3A + 2 * N3A_0)
            s = put_row_col_val(rowind, colind, -2.*value, row, row, val, s)

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
            s = put_row_col_val(rowind, colind, -2.*value, row, row-1, val, s)

            val = -non_kvsplit * N0E_0
            s = put_row_col_val(rowind, colind, -2.*value, row, row+2, val, s)

            val = + non_kvsplit * (N3E + 2 * N3E_0)
            s = put_row_col_val(rowind, colind, -2.*value, row, row, val, s)
        end

    end
end
