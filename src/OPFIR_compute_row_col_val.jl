function OPFIR_compute_row_col_val(rowind, colind, value, para, N::Array)
    s::Int64 = 1;
    row::Int64 = 0;
    col::Int64 = 0;
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib;
    unknown_layer::Int64 = layer_unknown;

    for vi::Int64 in 1:para.num_freq
        for ri::Int64 in 1:para.num_layers

            # oscil: ground; vib: J3
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 1;
            ## J3 -> J4, and thermal pool (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = -para.k12_G - para.k1a;
            s = s+1;
            ## J4 -> J3:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k21_G;
            s = s+1;

            # oscil: ground; vib: J4
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 2;
            ## J4 -> J5, J3 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k23_G - para.k21_G - para.k2a;
            s = s+1;
            ## J5, J3 -> J4:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k32_G;
            s = s+1;
            col = row -1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k12_G;
            s = s+1;

            ### pumping:
            col = row;
            rowind[s] = row;
            colind[s] = col;
            value[s] = - para.pumpR[vi];
            s = s+1;
            col = (ri-1)*unknown_layer + para.num_freq*n_rot + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = - para.pumpR[vi] * C4L * para.gauss_dist[vi];
            s = s+1;

            col = (ri-1)*unknown_layer + (vi-1)*n_rot + n_rot/2 + 3;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.pumpR[vi] * g_L/g_U;
            s = s+1;
            col = (ri-1)*unknown_layer + para.num_freq*n_rot + 2;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.pumpR[vi] * C5U * para.gauss_dist[vi] * g_L/g_U;
            s = s+1;

            # oscil: ground; vib: J5
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 3;
            ## J5 -> J6, J4 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k34_G - para.k32_G - para.k3a;
            s = s+1;
            ## J6, J4 -> J5:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k43_G;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k23_G;
            s = s+1;

            # oscil: ground; vib: J6
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 4;
            ## J6 -> J7, J5 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k45_G - para.k43_G - para.k4a;
            s = s+1;
            ## J7, J5 -> J6:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k54_G;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k34_G;
            s = s+1;

            # oscil: ground; vib: J7
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 5;
            ## J7 -> J8, J6 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k56_G - para.k54_G - para.k5a;
            s = s+1;
            ## J8, J6 -> J7:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k65_G;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k45_G;
            s = s+1;

            # oscil: ground; vib: J8
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 6;
            ## J8 -> J9, J7 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k67_G - para.k65_G - para.k6a;
            s = s+1;
            ## J9, J7 -> J8:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k76_G;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k56_G;
            s = s+1;

            # oscil: ground; vib: J9
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 7;
            ## J9 -> J10, J8 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k78_G - para.k76_G - para.k7a;
            s = s+1;
            ## J10, J8 -> J9:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k87_G;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k67_G;
            s = s+1;

            # oscil: ground; vib: J10
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 8;
            ## J10 -> J11, J9 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k89_G - para.k87_G - para.k8a;
            s = s+1;
            ## J11, J9 -> J11:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k98_G;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k78_G;
            s = s+1;

            # oscil: ground; vib: J11
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 9;
            ## J11 -> J10 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k98_G - para.k9a;
            s = s+1;
            ## J10 -> J11:
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k89_G;
            s = s+1;
        end
    end
    for vi::Int64 in 1:para.num_freq
        for ri::Int64 in 1:para.num_layers
            # oscil: 3; vib: J3
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 9 + 1;
            ## J3 -> J4 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = -para.k12_3 - para.k10a;
            s = s+1;
            ## J4 -> J3:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k21_3;
            s = s+1;

            # oscil: 3; vib: J4
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 9 + 2;
            ## J4 -> J5, J3 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k23_3 - para.k21_3 - para.k11a;
            s = s+1;
            ## J5, J3 -> J4:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k32_3;
            s = s+1;
            col = row -1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k12_3;
            s = s+1;

            # oscil: 3; vib: J5
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 9 + 3;
            ## J5 -> J6, J4 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k34_3 - para.k32_3 - para.k12a;
            s = s+1;
            ## J6, J4 -> J5:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k43_3;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k23_3;
            s = s+1;

            ### pumping:
            col = (ri-1)*unknown_layer + (vi-1)*n_rot + 2;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.pumpR[vi];
            s = s+1;
            col = (ri-1)*unknown_layer + para.num_freq*n_rot + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.pumpR[vi] * C4L * para.gauss_dist[vi];
            s = s+1;

            col = row;
            rowind[s] = row;
            colind[s] = col;
            value[s] = - para.pumpR[vi] * g_L/g_U;
            s = s+1;
            col = (ri-1)*unknown_layer + para.num_freq*n_rot + 2;
            rowind[s] = row;
            colind[s] = col;
            value[s] = - para.pumpR[vi] * C5L * para.gauss_dist[vi] * g_L/g_U;
            s = s+1;

            # oscil: 3; vib: J6
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 9 + 4;
            ## J6 -> J7, J5 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k45_3 - para.k43_3 - para.k13a;
            s = s+1;
            ## J7, J5 -> J6:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k54_3;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k34_3;
            s = s+1;

            # oscil: 3; vib: J7
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 9 + 5;
            ## J7 -> J8, J6 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k56_3 - para.k54_3 - para.k14a;
            s = s+1;
            ## J8, J6 -> J7:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k65_3;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k45_3;
            s = s+1;

            # oscil: 3; vib: J8
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 9 + 6;
            ## J8 -> J9, J7 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k67_3 - para.k65_3 - para.k15a;
            s = s+1;
            ## J9, J7 -> J8:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k76_3;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k56_3;
            s = s+1;

            # oscil: 3; vib: J9
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 9 + 7;
            ## J9 -> J10, J8 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k78_3 - para.k76_3 - para.k16a;
            s = s+1;
            ## J10, J8 -> J9:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k87_3;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k67_3;
            s = s+1;

            # oscil: 3; vib: J10
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 9 + 8;
            ## J10 -> J11, J9 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k89_3 - para.k87_3 - para.k17a;
            s = s+1;
            ## J11, J9 -> J11:
            col = row + 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k98_3;
            s = s+1;
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k78_3;
            s = s+1;

            # oscil: 3; vib: J11
            row = (ri-1)*unknown_layer + (vi-1)*n_rot + 9 + 9;
            ## J11 -> J10 (diagonal):
            rowind[s] = row;
            colind[s] = row;
            value[s] = - para.k98_3 - para.k18a;
            s = s+1;
            ## J10 -> J11:
            col = row - 1;
            rowind[s] = row;
            colind[s] = col;
            value[s] = para.k89_3;
            s = s+1;
        end
    end

    #### Thermal pools: ##########
    ### from rotational levels to thermal pools ###
    for ri in 1:para.num_layers
        # V0 A type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 1;
        for vi in 1:para.num_freq
            for k in 1:n_rot/2
                col = (ri-1)*unknown_layer + (vi-1)*n_rot + k;
                rowind[s] = row;
                colind[s] = col;
                value[s] = para.k1a;
                s = s+1;
            end
        end

        # V3 A type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 2;
        for vi in 1:para.num_freq
            for k in 1:9
                col = (ri-1)*unknown_layer + (vi-1)*n_rot + n_rot/2 + k;
                rowind[s] = row;
                colind[s] = col;
                value[s] = para.k1a;
                s = s+1;
            end
        end
    end

    ### transitions between thermal pools ###
    ## A type:
    for ri in 1:para.num_layers
        # V0 A type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 1;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.kwall[ri] + f_G * para.kwall[ri];
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_G * para.kwall[ri];
        s = s+1;

        col = row + 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_G * para.kwall[ri];
        s = s+1;

        col = row + 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_G * para.kwall[ri];
        s = s+1;

        col = row + 4;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_G * para.kwall[ri];
        s = s+1;

        col = row + 5;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_G * para.kwall[ri];
        s = s+1;

        # V3 A type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 2;
        rowind[s] = row;
        colind[s] = row;
        value[s] = - para.kwall[ri] + f_3 * para.kwall[ri];
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_3 * para.kwall[ri];
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_3 * para.kwall[ri];
        s = s+1;

        col = row + 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_3 * para.kwall[ri];
        s = s+1;

        col = row + 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_3 * para.kwall[ri];
        s = s+1;

        col = row + 4;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_3 * para.kwall[ri];
        s = s+1;

        #### transtion between V3 and V6 ####
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.k36;
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = +para.k63;
        s = s+1;

        # V6 A type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 3;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.kwall[ri] + f_6 * para.kwall[ri];
        s = s+1;

        col = row - 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_6 * para.kwall[ri];
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_6 * para.kwall[ri];
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_6 * para.kwall[ri];
        s = s+1;

        col = row + 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_6 * para.kwall[ri];
        s = s+1;

        col = row + 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_6 * para.kwall[ri];
        s = s+1;

        #### transtion between V3 and V6: ####
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.k63;
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = para.k36;
        s = s+1;

        # V23 A type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 4;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.kwall[ri] + f_23 * para.kwall[ri];
        s = s+1;

        col = row - 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_23 * para.kwall[ri];
        s = s+1;

        col = row - 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_23 * para.kwall[ri];
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_23 * para.kwall[ri];
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_23 * para.kwall[ri];
        s = s+1;

        col = row + 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_23 * para.kwall[ri];
        s = s+1;

        #### transitions between V23 and V36 ####
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.k2336;
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = para.k3623;
        s = s+1;

        # V36 A type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 5;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.kwall[ri] + f_36 * para.kwall[ri];
        s = s+1;

        col = row - 4;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_36 * para.kwall[ri];
        s = s+1;

        col = row - 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_36 * para.kwall[ri];
        s = s+1;

        col = row - 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_36 * para.kwall[ri];
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_36 * para.kwall[ri];
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_36 * para.kwall[ri];
        s = s+1;

        #### transition between V36 and V23, V26 ####
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.k3623 - para.k3626;
        s = s+1;

        rowind[s] = row;
        colind[s] = row - 1;
        value[s] = +para.k2336;
        s = s+1;

        rowind[s] = row;
        colind[s] = row + 1;
        value[s] = +para.k2636;
        s = s+1;

        # V26 A type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 6;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.kwall[ri] + f_26 * para.kwall[ri];
        s = s+1;

        col = row - 5;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_26 * para.kwall[ri];
        s = s+1;

        col = row - 4;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_26 * para.kwall[ri];
        s = s+1;

        col = row - 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_26 * para.kwall[ri];
        s = s+1;

        col = row - 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_26 * para.kwall[ri];
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_26 * para.kwall[ri];
        s = s+1;

        #### transition between V26 and V36: ####
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.k2636;
        s = s+1;

        rowind[s] = row;
        colind[s] = row - 1;
        value[s] = para.k3626;
        s = s+1;
    end

    ## E type:
    for ri in 1:para.num_layers
        # V0 E type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 7;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.kwall[ri] + f_G * para.kwall[ri];
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_G * para.kwall[ri];
        s = s+1;

        col = row + 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_G * para.kwall[ri];
        s = s+1;

        col = row + 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_G * para.kwall[ri];
        s = s+1;

        col = row + 4;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_G * para.kwall[ri];
        s = s+1;

        col = row + 5;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_G * para.kwall[ri];
        s = s+1;

        # V3 E type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 8;
        rowind[s] = row;
        colind[s] = row;
        value[s] = - para.kwall[ri] + f_3 * para.kwall[ri];
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_3 * para.kwall[ri];
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_3 * para.kwall[ri];
        s = s+1;

        col = row + 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_3 * para.kwall[ri];
        s = s+1;

        col = row + 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_3 * para.kwall[ri];
        s = s+1;

        col = row + 4;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_3 * para.kwall[ri];
        s = s+1;

        #### transtion between V3 and V6 ####
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.k36;
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = +para.k63;
        s = s+1;

        # V6 E type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 9;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.kwall[ri] + f_6 * para.kwall[ri];
        s = s+1;

        col = row - 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_6 * para.kwall[ri];
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_6 * para.kwall[ri];
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_6 * para.kwall[ri];
        s = s+1;

        col = row + 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_6 * para.kwall[ri];
        s = s+1;

        col = row + 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_6 * para.kwall[ri];
        s = s+1;

        #### transtion between V3 and V6: ####
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.k63;
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = para.k36;
        s = s+1;

        # V23 E type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 10;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.kwall[ri] + f_23 * para.kwall[ri];
        s = s+1;

        col = row - 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_23 * para.kwall[ri];
        s = s+1;

        col = row - 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_23 * para.kwall[ri];
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_23 * para.kwall[ri];
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_23 * para.kwall[ri];
        s = s+1;

        col = row + 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_23 * para.kwall[ri];
        s = s+1;

        #### transitions between V23 and V36 ####
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.k2336;
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = para.k3623;
        s = s+1;

        # V36 E type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 11;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.kwall[ri] + f_36 * para.kwall[ri];
        s = s+1;

        col = row - 4;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_36 * para.kwall[ri];
        s = s+1;

        col = row - 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_36 * para.kwall[ri];
        s = s+1;

        col = row - 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_36 * para.kwall[ri];
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_36 * para.kwall[ri];
        s = s+1;

        col = row + 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_36 * para.kwall[ri];
        s = s+1;

        #### transition between V36 and V23, V26 ####
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.k3623 - para.k3626;
        s = s+1;

        rowind[s] = row;
        colind[s] = row - 1;
        value[s] = +para.k2336;
        s = s+1;

        rowind[s] = row;
        colind[s] = row + 1;
        value[s] = +para.k2636;
        s = s+1;

        # V26 E type:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 12;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.kwall[ri] + f_26 * para.kwall[ri];
        s = s+1;

        col = row - 5;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_26 * para.kwall[ri];
        s = s+1;

        col = row - 4;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_26 * para.kwall[ri];
        s = s+1;

        col = row - 3;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_26 * para.kwall[ri];
        s = s+1;

        col = row - 2;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_26 * para.kwall[ri];
        s = s+1;

        col = row - 1;
        rowind[s] = row;
        colind[s] = col;
        value[s] = f_26 * para.kwall[ri];
        s = s+1;

        #### transition between V26 and V36: ####
        rowind[s] = row;
        colind[s] = row;
        value[s] = -para.k2636;
        s = s+1;

        rowind[s] = row;
        colind[s] = row - 1;
        value[s] = para.k3626;
        s = s+1;
    end

    ########################## add the diffusion term ####################################

    for ri::Int64 in 2:para.num_layers-1
        index_diffu = layer_unknown * (ri-1) + 1 : layer_unknown * ri;
        ## close diffusion for the 2 E type levels:
        # index_diffu = layer_unknown * (ri-1) + 1 : layer_unknown * ri - 2;

        for k in index_diffu
            rowind[s] = k;
            colind[s] = k - layer_unknown;
            value[s] = para.D*(-1.0/(2*para.Δr*para.r_int[ri]) + 1.0/(para.Δr)^2);
            s = s+1;

            rowind[s] = k;
            colind[s] = k;
            value[s] = para.D*(-2.0/(para.Δr)^2);
            s = s+1;

            rowind[s] = k;
            colind[s] = k + layer_unknown;
            value[s] = para.D*(1.0/(2*para.Δr*para.r_int[ri]) + 1.0/(para.Δr)^2);
            s = s+1;
        end
    end

    for k in 1:layer_unknown
    ## close the diffusion for the 2 E type levels
    # for k in 1:layer_unknown - 2
        rowind[s] = k;
        colind[s] = k;
        value[s] = para.D*(-1.0/(para.Δr)^2 - 1.0/2.0/para.Δr/para.r_int[1]);
        s = s+1;

        rowind[s] = k;
        colind[s] = k + layer_unknown;
        value[s] = para.D*(1.0/(para.Δr)^2 + 1.0/2/para.Δr/para.r_int[1]);
        s = s+1;
    end

    ################ question: BC for the wall? #######
    for k in layer_unknown*(para.num_layers-1)+1:layer_unknown*para.num_layers
    ## close diffusion for the 2 E type levels:
    # for k in layer_unknown*(para.num_layers-1)+1:layer_unknown*para.num_layers-2
        rowind[s] = k;
        colind[s] = k - layer_unknown;
        value[s] = para.D*(-1.0/(2*para.Δr*para.r_int[end]) + 1.0/(para.Δr)^2);
        s = s+1;

        rowind[s] = k;
        colind[s] = k;
        value[s] = para.D*(-2.0/(para.Δr)^2) +
                   para.D*(1.0/(2*para.Δr*para.r_int[end]) + 1.0/(para.Δr)^2);
        s = s+1;
    end

    ############## add the V-swap process V0A + V3E <-> V_3A + V0E> ##########
    N0A_0 = para.ntotal * f_G/2;
    N0E_0 = N0A_0;
    N3A_0 = para.ntotal * f_3/2;
    N3E_0 = N3A_0;
    for ri::Int64 in 1:para.num_layers

        offset = (ri-1)*unknown_layer + para.num_freq*n_rot;
        N0A = N[offset + 1];
        N3A = N[offset + 2];
        N0E = N[offset + 6 + 1];
        N3E = N[offset + 6 + 2];
        ## V0A:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 1;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -kvs * N3E - kvs * N3E_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row + 6 + 1; # V3E Level
        value[s] = -kvs * N0A_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row + 6; # based on N0E;
        value[s] = + kvs * N3A + kvs * N3A_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row + 1;
        value[s] = + kvs * N0E_0;
        s = s+1;

        ## V3E
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 6 + 2;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -kvs * N0A - kvs * N0A_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row - 6 - 1;
        value[s] = -kvs * N3E_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row - 1;
        value[s] = + kvs * N3A_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row - 6;
        value[s] = + kvs * N0E + kvs * N0E_0;
        s = s+1;

        ## V0E:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 6 + 1;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -kvs * N3A - kvs * N3A_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row - 6 + 1;
        value[s] = -kvs * N0E_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row + 1;
        value[s] = +kvs * N0A_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row - 6;
        value[s] = + kvs * N3E + kvs * N3E_0;
        s = s+1;

        ## V3A:
        row = (ri-1)*unknown_layer + para.num_freq*n_rot + 2;
        rowind[s] = row;
        colind[s] = row;
        value[s] = -kvs * N0E - kvs * N0E_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row + 6 - 1;
        value[s] = -kvs * N3A_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row + 6;
        value[s] = + kvs * N0A + kvs * N0A_0;
        s = s+1;

        rowind[s] = row;
        colind[s] = row - 1;
        value[s] = + kvs * N3E_0;
        s = s+1;
    end

end
