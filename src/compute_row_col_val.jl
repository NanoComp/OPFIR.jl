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
            val = -p.k12_G - p.k1a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J4 -> J3:
            val = p.k21_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # oscil: ground vib: J4
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 2
            ## J4 -> J5, J3 (diagonal):
            val = - p.k23_G - p.k21_G - p.k2a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J5, J3 -> J4:
            val = p.k32_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k12_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
            ### pumping:
            val = -p.pumpR[vi]
            s = put_row_col_val(rowind, colind, value, row, row, val, s)

            val = - p.pumpR[vi] * p.C4L * p.gauss_dist[vi]
            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + 3
            val = p.pumpR[vi] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
            val = p.pumpR[vi] * p.C5U * p.gauss_dist[vi] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            # oscil: ground vib: J5
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 3
            ## J5 -> J6, J4 (diagonal):
            val = - p.k34_G - p.k32_G - p.k3a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J6, J4 -> J5:
            val = p.k43_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k23_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J6
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 4
            ## J6 -> J7, J5 (diagonal):
            val = - p.k45_G - p.k43_G - p.k4a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J7, J5 -> J6:
            val = p.k54_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k34_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J7
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 5
            ## J7 -> J8, J6 (diagonal):
            val = - p.k56_G - p.k54_G - p.k5a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J8, J6 -> J7:
            val = p.k65_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k45_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J8
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 6
            ## J8 -> J9, J7 (diagonal):
            val = - p.k67_G - p.k65_G - p.k6a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J9, J7 -> J8:
            val = p.k76_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k56_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J9
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 7
            ## J9 -> J10, J8 (diagonal):
            val = - p.k78_G - p.k76_G - p.k7a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J10, J8 -> J9:
            val = p.k87_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k67_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J10
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 8
            ## J10 -> J11, J9 (diagonal):
            val = - p.k89_G - p.k87_G - p.k8a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J11, J9 -> J11:
            val = p.k98_G
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k78_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: ground vib: J11
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9
            ## J11 -> J10 (diagonal):
            val = - p.k98_G - p.k9a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J10 -> J11:
            val = p.k89_G
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        end
    end
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
            # oscil: 3 vib: J3
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 1
            ## J3 -> J4 (diagonal):
            val = -p.k12_3 - p.k10a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J4 -> J3:
            val = p.k21_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

            # oscil: 3 vib: J4
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 2
            ## J4 -> J5, J3 (diagonal):
            val = - p.k23_3 - p.k21_3 - p.k11a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J5, J3 -> J4:
            val = p.k32_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k12_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J5
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 3
            ## J5 -> J6, J4 (diagonal):
            val = - p.k34_3 - p.k32_3 - p.k12a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J6, J4 -> J5:
            val = p.k43_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k23_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            ### pumping:
            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 2
            val = p.pumpR[vi]
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
            val = p.pumpR[vi] * p.C4L * p.gauss_dist[vi]
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = row
            val = - p.pumpR[vi] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
            val = - p.pumpR[vi] * p.C5L * p.gauss_dist[vi] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            # oscil: 3 vib: J6
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 4
            ## J6 -> J7, J5 (diagonal):
            val = - p.k45_3 - p.k43_3 - p.k13a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J7, J5 -> J6:
            val = p.k54_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k34_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J7
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 5
            ## J7 -> J8, J6 (diagonal):
            val = - p.k56_3 - p.k54_3 - p.k14a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J8, J6 -> J7:
            val = p.k65_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k45_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J8
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 6
            ## J8 -> J9, J7 (diagonal):
            val = - p.k67_3 - p.k65_3 - p.k15a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J9, J7 -> J8:
            val = p.k76_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k56_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J9
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 7
            ## J9 -> J10, J8 (diagonal):
            val = - p.k78_3 - p.k76_3 - p.k16a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J10, J8 -> J9:
            val = p.k87_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k67_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J10
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 8
            ## J10 -> J11, J9 (diagonal):
            val = - p.k89_3 - p.k87_3 - p.k17a
            s = put_row_col_val(rowind, colind, value, row, row, val, s)
            ## J11, J9 -> J11:
            val = p.k98_3
            s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
            val = p.k78_3
            s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

            # oscil: 3 vib: J11
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 9 + 9
            ## J11 -> J10 (diagonal):
            val = - p.k98_3 - p.k18a
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
                val = p.k1a
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
        end

        # V3 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
        for vi in 1:p.num_freq
            for k in 1:9
                col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + k
                val = p.k1a
                s = put_row_col_val(rowind, colind, value, row, col, val, s)
            end
        end
    end

    ### transitions between thermal pools ###
    ## A type:
    for ri in 1:p.num_layers
        # V0 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1
        val =  -p.kwall[ri] + p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        val = p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+2, val, s)
        val = p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+3, val, s)
        val = p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+4, val, s)
        val = p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+5, val, s)

        # V3 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
        val = - p.kwall[ri] + p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        val = p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        val = p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+2, val, s)
        val = p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+3, val, s)
        val = p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+4, val, s)

        #### transtion between V3 and V6 ####
        val = -p.k36
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = +p.k63
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        # V6 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 3
        val = -p.kwall[ri] + p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-2, val, s)
        val = p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        val = p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        val = p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+2, val, s)
        val = p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+3, val, s)

        #### transtion between V3 and V6: ####
        val = -p.k63
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = +p.k36
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

        # V23 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 4
        val = -p.kwall[ri] + p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-3, val, s)
        val = p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-2, val, s)
        val = p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        val = p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        val = p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+2, val, s)

        #### transitions between V23 and V36 ####
        val = -p.k2336
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = +p.k3623
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        # V36 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 5
        val = -p.kwall[ri] + p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-4, val, s)
        val = p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-3, val, s)
        val = p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-2, val, s)
        val = p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        val = p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        #### transition between V36 and V23, V26 ####
        val = -p.k3623 - p.k3626
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = +p.k2336
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        val = +p.k2636
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        # V26 A type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 6
        val = -p.kwall[ri] + p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-5, val, s)
        val = p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-4, val, s)
        val = p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-3, val, s)
        val = p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-2, val, s)
        val = p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

        #### transition between V26 and V36: ####
        val = -p.k2636
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = +p.k3626
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
    end

    ## E type:
    for ri in 1:p.num_layers
        # V0 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 7
        val = -p.kwall[ri] + p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        val = p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+2, val, s)
        val = p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+3, val, s)
        val = p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+4, val, s)
        val = p.f_G * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+5, val, s)

        # V3 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 8
        val = - p.kwall[ri] + p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        val = p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        val = p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+2, val, s)
        val = p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+3, val, s)
        val = p.f_3 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+4, val, s)

        #### transtion between V3 and V6 ####
        val = -p.k36
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = +p.k63
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        # V6 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 9
        val = -p.kwall[ri] + p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-2, val, s)
        val = p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        val = p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        val = p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+2, val, s)
        val = p.f_6 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+3, val, s)

        #### transtion between V3 and V6: ####
        val = -p.k63
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = +p.k36
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

        # V23 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 10
        val = -p.kwall[ri] + p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-3, val, s)
        val = p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-2, val, s)
        val = p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        val = p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)
        val = p.f_23 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+2, val, s)

        #### transitions between V23 and V36 ####
        val = -p.k2336
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = +p.k3623
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        # V36 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 11
        val = -p.kwall[ri] + p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-4, val, s)
        val = p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-3, val, s)
        val = p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-2, val, s)
        val = p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        val = p.f_36 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        #### transition between V36 and V23, V26 ####
        val = -p.k3623 - p.k3626
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = +p.k2336
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
        val = +p.k2636
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        # V26 E type:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 12
        val = -p.kwall[ri] + p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-5, val, s)
        val = p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-4, val, s)
        val = p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-3, val, s)
        val = p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-2, val, s)
        val = p.f_26 * p.kwall[ri]
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

        #### transition between V26 and V36: ####
        val = -p.k2636
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
        val = +p.k3626
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
    end

    ########################## add the diffusion term ####################################

    for ri in 2:p.num_layers-1
        index_diffu = p.layer_unknown * (ri-1) + 1 : p.layer_unknown * ri
        ## close diffusion for the 2 E type levels:
        # index_diffu = p.layer_unknown * (ri-1) + 1 : p.layer_unknown * ri - 2

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

    for k in 1:p.layer_unknown
    ## close the diffusion for the 2 E type levels
    # for k in 1:p.layer_unknown - 2
        row = k
        val = p.D*(-1.0/(p.Δr)^2 - 1.0/2.0/p.Δr/p.r_int[1])
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        col = row + p.layer_unknown
        val = p.D*(1.0/(p.Δr)^2 + 1.0/2/p.Δr/p.r_int[1])
        s = put_row_col_val(rowind, colind, value, row, col, val, s)
    end

    ################ question: BC for the wall? #######
    for k in p.layer_unknown*(p.num_layers-1)+1:p.layer_unknown*p.num_layers
    ## close diffusion for the 2 E type levels:
    # for k in p.layer_unknown*(p.num_layers-1)+1:p.layer_unknown*p.num_layers-2
        row = k
        col = row - p.layer_unknown
        val = p.D*(-1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
        s = put_row_col_val(rowind, colind, value, row, col, val, s)

        val = p.D*(-2.0/(p.Δr)^2) + p.D*(1.0/(2*p.Δr*p.r_int[end]) + 1.0/(p.Δr)^2)
        s = put_row_col_val(rowind, colind, value, row, row, val, s)
    end

    ############## add the V-swap process V0A + V3E <-> V_3A + V0E> ##########
    N0A_0 = p.ntotal * p.f_G/2
    N0E_0 = N0A_0
    N3A_0 = p.ntotal * p.f_3/2
    N3E_0 = N3A_0
    for ri in 1:p.num_layers
        offset = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot
        N0A = sol_0[offset + 1]
        N3A = sol_0[offset + 2]
        N0E = sol_0[offset + 6 + 1]
        N3E = sol_0[offset + 6 + 2]
        ## V0A:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 1

        val = -p.kvs * N3E - p.kvs * N3E_0
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = -p.kvs * N0A_0
        s = put_row_col_val(rowind, colind, value, row, row+6+1, val, s)

        val = + p.kvs * N3A + p.kvs * N3A_0
        s = put_row_col_val(rowind, colind, value, row, row+6, val, s)

        val = + p.kvs * N0E_0
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        ## V3E
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 6 + 2

        val = -p.kvs * N0A - p.kvs * N0A_0
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = -p.kvs * N3E_0
        s = put_row_col_val(rowind, colind, value, row, row-6-1, val, s)

        val = + p.kvs * N3A_0
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)

        val = + p.kvs * N0E + p.kvs * N0E_0
        s = put_row_col_val(rowind, colind, value, row, row-6, val, s)

        ## V0E:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 6 + 1
        val = -p.kvs * N3A - p.kvs * N3A_0
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = -p.kvs * N0E_0
        s = put_row_col_val(rowind, colind, value, row, row-6+1, val, s)

        val = +p.kvs * N0A_0
        s = put_row_col_val(rowind, colind, value, row, row+1, val, s)

        val = + p.kvs * N3E + p.kvs * N3E_0
        s = put_row_col_val(rowind, colind, value, row, row-6, val, s)

        ## V3A:
        row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2

        val = -p.kvs * N0E - p.kvs * N0E_0
        s = put_row_col_val(rowind, colind, value, row, row, val, s)

        val = -p.kvs * N3A_0
        s = put_row_col_val(rowind, colind, value, row, row+6-1, val, s)

        val = + p.kvs * N0A + p.kvs * N0A_0
        s = put_row_col_val(rowind, colind, value, row, row+6, val, s)

        val = + p.kvs * N3E_0
        s = put_row_col_val(rowind, colind, value, row, row-1, val, s)
    end
end