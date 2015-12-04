function put_row_col_val(rowind, colind, value, row, col, val, s)
    rowind[s] = row
    colind[s] = col
    value[s] = val
    return (s+1)
end

function compute_row_col_val(rowind, colind, value, p, sol_0)
    s = 1
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
            # oscil: ground vib: J3
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 1
            ## J3 -> J4, and thermal pool (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = -p.k12_G - p.k1a
            s = s+1
            ## J4 -> J3:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k21_G
            s = s+1

            # oscil: ground vib: J4
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 2
            ## J4 -> J5, J3 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k23_G - p.k21_G - p.k2a
            s = s+1
            ## J5, J3 -> J4:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k32_G
            s = s+1
            col = row -1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k12_G
            s = s+1

            ### pumping:
            col = row
            rowind[s] = row
            colind[s] = col
            value[s] = - p.pumpR[vi]
            s = s+1
            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
            rowind[s] = row
            colind[s] = col
            value[s] = - p.pumpR[vi] * p.C4L * p.gauss_dist[vi]
            s = s+1

            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + 3
            rowind[s] = row
            colind[s] = col
            value[s] = p.pumpR[vi] * p.g_L/p.g_U
            s = s+1
            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
            rowind[s] = row
            colind[s] = col
            value[s] = p.pumpR[vi] * p.C5U * p.gauss_dist[vi] * p.g_L/p.g_U
            s = s+1

            # oscil: ground vib: J5
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 3
            ## J5 -> J6, J4 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k34_G - p.k32_G - p.k3a
            s = s+1
            ## J6, J4 -> J5:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k43_G
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k23_G
            s = s+1

            # oscil: ground vib: J6
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 4
            ## J6 -> J7, J5 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k45_G - p.k43_G - p.k4a
            s = s+1
            ## J7, J5 -> J6:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k54_G
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k34_G
            s = s+1

            # oscil: ground vib: J7
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 5
            ## J7 -> J8, J6 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k56_G - p.k54_G - p.k5a
            s = s+1
            ## J8, J6 -> J7:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k65_G
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k45_G
            s = s+1

            # oscil: ground vib: J8
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 6
            ## J8 -> J9, J7 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k67_G - p.k65_G - p.k6a
            s = s+1
            ## J9, J7 -> J8:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k76_G
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k56_G
            s = s+1

            # oscil: ground vib: J9
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 7
            ## J9 -> J10, J8 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k78_G - p.k76_G - p.k7a
            s = s+1
            ## J10, J8 -> J9:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k87_G
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k67_G
            s = s+1

            # oscil: ground vib: J10
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 8
            ## J10 -> J11, J9 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k89_G - p.k87_G - p.k8a
            s = s+1
            ## J11, J9 -> J11:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k98_G
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k78_G
            s = s+1

            # oscil: ground vib: J11
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9
            ## J11 -> J10 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k98_G - p.k9a
            s = s+1
            ## J10 -> J11:
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k89_G
            s = s+1
        end
    end
    for vi::Int64 in 1:p.num_freq
        for ri::Int64 in 1:p.num_layers
            # oscil: 3 vib: J3
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 1
            ## J3 -> J4 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = -p.k12_3 - p.k10a
            s = s+1
            ## J4 -> J3:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k21_3
            s = s+1

            # oscil: 3 vib: J4
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 2
            ## J4 -> J5, J3 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k23_3 - p.k21_3 - p.k11a
            s = s+1
            ## J5, J3 -> J4:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k32_3
            s = s+1
            col = row -1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k12_3
            s = s+1

            # oscil: 3 vib: J5
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 3
            ## J5 -> J6, J4 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k34_3 - p.k32_3 - p.k12a
            s = s+1
            ## J6, J4 -> J5:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k43_3
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k23_3
            s = s+1

            ### pumping:
            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 2
            rowind[s] = row
            colind[s] = col
            value[s] = p.pumpR[vi]
            s = s+1
            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.pumpR[vi] * p.C4L * p.gauss_dist[vi]
            s = s+1

            col = row
            rowind[s] = row
            colind[s] = col
            value[s] = - p.pumpR[vi] * p.g_L/p.g_U
            s = s+1
            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
            rowind[s] = row
            colind[s] = col
            value[s] = - p.pumpR[vi] * p.C5L * p.gauss_dist[vi] * p.g_L/p.g_U
            s = s+1

            # oscil: 3 vib: J6
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 4
            ## J6 -> J7, J5 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k45_3 - p.k43_3 - p.k13a
            s = s+1
            ## J7, J5 -> J6:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k54_3
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k34_3
            s = s+1

            # oscil: 3 vib: J7
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 5
            ## J7 -> J8, J6 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k56_3 - p.k54_3 - p.k14a
            s = s+1
            ## J8, J6 -> J7:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k65_3
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k45_3
            s = s+1

            # oscil: 3 vib: J8
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 6
            ## J8 -> J9, J7 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k67_3 - p.k65_3 - p.k15a
            s = s+1
            ## J9, J7 -> J8:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k76_3
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k56_3
            s = s+1

            # oscil: 3 vib: J9
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 7
            ## J9 -> J10, J8 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k78_3 - p.k76_3 - p.k16a
            s = s+1
            ## J10, J8 -> J9:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k87_3
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k67_3
            s = s+1

            # oscil: 3 vib: J10
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 8
            ## J10 -> J11, J9 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k89_3 - p.k87_3 - p.k17a
            s = s+1
            ## J11, J9 -> J11:
            col = row + 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k98_3
            s = s+1
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k78_3
            s = s+1

            # oscil: 3 vib: J11
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 9
            ## J11 -> J10 (diagonal):
            rowind[s] = row
            colind[s] = row
            value[s] = - p.k98_3 - p.k18a
            s = s+1
            ## J10 -> J11:
            col = row - 1
            rowind[s] = row
            colind[s] = col
            value[s] = p.k89_3
            s = s+1
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
                rowind[s] = row
                colind[s] = col
                value[s] = p.k1a
                s = s+1
            end
        end

        # V3 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
        for vi in 1:p.num_freq
            for k in 1:9
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + k
                rowind[s] = row
                colind[s] = col
                value[s] = p.k1a
                s = s+1
            end
        end
    end

    ### transitions between thermal pools ###
    ## A type:
    for ri in 1:p.num_layers
        # V0 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kwall[ri] + p.f_G * p.kwall[ri]
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_G * p.kwall[ri]
        s = s+1

        col = row + 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_G * p.kwall[ri]
        s = s+1

        col = row + 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_G * p.kwall[ri]
        s = s+1

        col = row + 4
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_G * p.kwall[ri]
        s = s+1

        col = row + 5
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_G * p.kwall[ri]
        s = s+1

        # V3 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
        rowind[s] = row
        colind[s] = row
        value[s] = - p.kwall[ri] + p.f_3 * p.kwall[ri]
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_3 * p.kwall[ri]
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_3 * p.kwall[ri]
        s = s+1

        col = row + 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_3 * p.kwall[ri]
        s = s+1

        col = row + 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_3 * p.kwall[ri]
        s = s+1

        col = row + 4
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_3 * p.kwall[ri]
        s = s+1

        #### transtion between V3 and V6 ####
        rowind[s] = row
        colind[s] = row
        value[s] = -p.k36
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = +p.k63
        s = s+1

        # V6 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 3
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kwall[ri] + p.f_6 * p.kwall[ri]
        s = s+1

        col = row - 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_6 * p.kwall[ri]
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_6 * p.kwall[ri]
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_6 * p.kwall[ri]
        s = s+1

        col = row + 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_6 * p.kwall[ri]
        s = s+1

        col = row + 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_6 * p.kwall[ri]
        s = s+1

        #### transtion between V3 and V6: ####
        rowind[s] = row
        colind[s] = row
        value[s] = -p.k63
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.k36
        s = s+1

        # V23 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 4
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kwall[ri] + p.f_23 * p.kwall[ri]
        s = s+1

        col = row - 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_23 * p.kwall[ri]
        s = s+1

        col = row - 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_23 * p.kwall[ri]
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_23 * p.kwall[ri]
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_23 * p.kwall[ri]
        s = s+1

        col = row + 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_23 * p.kwall[ri]
        s = s+1

        #### transitions between V23 and V36 ####
        rowind[s] = row
        colind[s] = row
        value[s] = -p.k2336
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.k3623
        s = s+1

        # V36 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 5
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kwall[ri] + p.f_36 * p.kwall[ri]
        s = s+1

        col = row - 4
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_36 * p.kwall[ri]
        s = s+1

        col = row - 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_36 * p.kwall[ri]
        s = s+1

        col = row - 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_36 * p.kwall[ri]
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_36 * p.kwall[ri]
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_36 * p.kwall[ri]
        s = s+1

        #### transition between V36 and V23, V26 ####
        rowind[s] = row
        colind[s] = row
        value[s] = -p.k3623 - p.k3626
        s = s+1

        rowind[s] = row
        colind[s] = row - 1
        value[s] = +p.k2336
        s = s+1

        rowind[s] = row
        colind[s] = row + 1
        value[s] = +p.k2636
        s = s+1

        # V26 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 6
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kwall[ri] + p.f_26 * p.kwall[ri]
        s = s+1

        col = row - 5
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_26 * p.kwall[ri]
        s = s+1

        col = row - 4
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_26 * p.kwall[ri]
        s = s+1

        col = row - 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_26 * p.kwall[ri]
        s = s+1

        col = row - 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_26 * p.kwall[ri]
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_26 * p.kwall[ri]
        s = s+1

        #### transition between V26 and V36: ####
        rowind[s] = row
        colind[s] = row
        value[s] = -p.k2636
        s = s+1

        rowind[s] = row
        colind[s] = row - 1
        value[s] = p.k3626
        s = s+1
    end

    ## E type:
    for ri in 1:p.num_layers
        # V0 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 7
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kwall[ri] + p.f_G * p.kwall[ri]
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_G * p.kwall[ri]
        s = s+1

        col = row + 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_G * p.kwall[ri]
        s = s+1

        col = row + 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_G * p.kwall[ri]
        s = s+1

        col = row + 4
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_G * p.kwall[ri]
        s = s+1

        col = row + 5
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_G * p.kwall[ri]
        s = s+1

        # V3 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 8
        rowind[s] = row
        colind[s] = row
        value[s] = - p.kwall[ri] + p.f_3 * p.kwall[ri]
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_3 * p.kwall[ri]
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_3 * p.kwall[ri]
        s = s+1

        col = row + 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_3 * p.kwall[ri]
        s = s+1

        col = row + 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_3 * p.kwall[ri]
        s = s+1

        col = row + 4
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_3 * p.kwall[ri]
        s = s+1

        #### transtion between V3 and V6 ####
        rowind[s] = row
        colind[s] = row
        value[s] = -p.k36
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = +p.k63
        s = s+1

        # V6 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 9
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kwall[ri] + p.f_6 * p.kwall[ri]
        s = s+1

        col = row - 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_6 * p.kwall[ri]
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_6 * p.kwall[ri]
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_6 * p.kwall[ri]
        s = s+1

        col = row + 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_6 * p.kwall[ri]
        s = s+1

        col = row + 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_6 * p.kwall[ri]
        s = s+1

        #### transtion between V3 and V6: ####
        rowind[s] = row
        colind[s] = row
        value[s] = -p.k63
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.k36
        s = s+1

        # V23 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 10
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kwall[ri] + p.f_23 * p.kwall[ri]
        s = s+1

        col = row - 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_23 * p.kwall[ri]
        s = s+1

        col = row - 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_23 * p.kwall[ri]
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_23 * p.kwall[ri]
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_23 * p.kwall[ri]
        s = s+1

        col = row + 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_23 * p.kwall[ri]
        s = s+1

        #### transitions between V23 and V36 ####
        rowind[s] = row
        colind[s] = row
        value[s] = -p.k2336
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.k3623
        s = s+1

        # V36 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 11
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kwall[ri] + p.f_36 * p.kwall[ri]
        s = s+1

        col = row - 4
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_36 * p.kwall[ri]
        s = s+1

        col = row - 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_36 * p.kwall[ri]
        s = s+1

        col = row - 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_36 * p.kwall[ri]
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_36 * p.kwall[ri]
        s = s+1

        col = row + 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_36 * p.kwall[ri]
        s = s+1

        #### transition between V36 and V23, V26 ####
        rowind[s] = row
        colind[s] = row
        value[s] = -p.k3623 - p.k3626
        s = s+1

        rowind[s] = row
        colind[s] = row - 1
        value[s] = +p.k2336
        s = s+1

        rowind[s] = row
        colind[s] = row + 1
        value[s] = +p.k2636
        s = s+1

        # V26 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 12
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kwall[ri] + p.f_26 * p.kwall[ri]
        s = s+1

        col = row - 5
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_26 * p.kwall[ri]
        s = s+1

        col = row - 4
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_26 * p.kwall[ri]
        s = s+1

        col = row - 3
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_26 * p.kwall[ri]
        s = s+1

        col = row - 2
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_26 * p.kwall[ri]
        s = s+1

        col = row - 1
        rowind[s] = row
        colind[s] = col
        value[s] = p.f_26 * p.kwall[ri]
        s = s+1

        #### transition between V26 and V36: ####
        rowind[s] = row
        colind[s] = row
        value[s] = -p.k2636
        s = s+1

        rowind[s] = row
        colind[s] = row - 1
        value[s] = p.k3626
        s = s+1
    end

    ########################## add the diffusion term ####################################

    for ri::Int64 in 2:p.num_layers-1
        index_diffu = p.layer_unknown * (ri-1) + 1 : p.layer_unknown * ri
        ## close diffusion for the 2 E type levels:
        # index_diffu = p.layer_unknown * (ri-1) + 1 : p.layer_unknown * ri - 2

        for k in index_diffu
            rowind[s] = k
            colind[s] = k - p.layer_unknown
            value[s] = p.D*(-1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
            s = s+1

            rowind[s] = k
            colind[s] = k
            value[s] = p.D*(-2.0/(p.Δr)^2)
            s = s+1

            rowind[s] = k
            colind[s] = k + p.layer_unknown
            value[s] = p.D*(1.0/(2*p.Δr*p.r_int[ri]) + 1.0/(p.Δr)^2)
            s = s+1
        end
    end

    for k in 1:p.layer_unknown
    ## close the diffusion for the 2 E type levels
    # for k in 1:p.layer_unknown - 2
        rowind[s] = k
        colind[s] = k
        value[s] = p.D*(-1.0/(p.Δr)^2 - 1.0/2.0/p.Δr/p.r_int[1])
        s = s+1

        rowind[s] = k
        colind[s] = k + p.layer_unknown
        value[s] = p.D*(1.0/(p.Δr)^2 + 1.0/2/p.Δr/p.r_int[1])
        s = s+1
    end

    ################ question: BC for the wall? #######
    for k in p.layer_unknown*(p.num_layers-1)+1:p.layer_unknown*p.num_layers
    ## close diffusion for the 2 E type levels:
    # for k in p.layer_unknown*(p.num_layers-1)+1:p.layer_unknown*p.num_layers-2
        rowind[s] = k
        colind[s] = k - p.layer_unknown
        value[s] = p.D*(-1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
        s = s+1

        rowind[s] = k
        colind[s] = k
        value[s] = p.D*(-2.0/(p.Δr)^2) +
                   p.D*(1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
        s = s+1
    end

    ############## add the V-swap process V0A + V3E <-> V_3A + V0E> ##########
    N0A_0 = p.ntotal * p.f_G/2
    N0E_0 = N0A_0
    N3A_0 = p.ntotal * p.f_3/2
    N3E_0 = N3A_0
    for ri::Int64 in 1:p.num_layers

        offset = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot
        N0A = sol_0[offset + 1]
        N3A = sol_0[offset + 2]
        N0E = sol_0[offset + 6 + 1]
        N3E = sol_0[offset + 6 + 2]
        ## V0A:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kvs * N3E - p.kvs * N3E_0
        s = s+1

        rowind[s] = row
        colind[s] = row + 6 + 1 # V3E Level
        value[s] = -p.kvs * N0A_0
        s = s+1

        rowind[s] = row
        colind[s] = row + 6 # based on N0E
        value[s] = + p.kvs * N3A + p.kvs * N3A_0
        s = s+1

        rowind[s] = row
        colind[s] = row + 1
        value[s] = + p.kvs * N0E_0
        s = s+1

        ## V3E
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 6 + 2
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kvs * N0A - p.kvs * N0A_0
        s = s+1

        rowind[s] = row
        colind[s] = row - 6 - 1
        value[s] = -p.kvs * N3E_0
        s = s+1

        rowind[s] = row
        colind[s] = row - 1
        value[s] = + p.kvs * N3A_0
        s = s+1

        rowind[s] = row
        colind[s] = row - 6
        value[s] = + p.kvs * N0E + p.kvs * N0E_0
        s = s+1

        ## V0E:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 6 + 1
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kvs * N3A - p.kvs * N3A_0
        s = s+1

        rowind[s] = row
        colind[s] = row - 6 + 1
        value[s] = -p.kvs * N0E_0
        s = s+1

        rowind[s] = row
        colind[s] = row + 1
        value[s] = +p.kvs * N0A_0
        s = s+1

        rowind[s] = row
        colind[s] = row - 6
        value[s] = + p.kvs * N3E + p.kvs * N3E_0
        s = s+1

        ## V3A:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
        rowind[s] = row
        colind[s] = row
        value[s] = -p.kvs * N0E - p.kvs * N0E_0
        s = s+1

        rowind[s] = row
        colind[s] = row + 6 - 1
        value[s] = -p.kvs * N3A_0
        s = s+1

        rowind[s] = row
        colind[s] = row + 6
        value[s] = + p.kvs * N0A + p.kvs * N0A_0
        s = s+1

        rowind[s] = row
        colind[s] = row - 1
        value[s] = + p.kvs * N3E_0
        s = s+1
    end

end
