function put_row_col_val(rowind, colind, value, row, col, val, s)
    rowind[s] = row
    colind[s] = col
    value[s] = val
    return (s+1)
end

function compute_row_col_val(rowind, colind, value, p, sol_0)
    s = 1
    A = 1000
    B = 1000.
    # rotational levels in V0 and V3
    for vi in 1:p.num_freq
        for ri in 1:p.num_layers
            for li in 1:p.n_rot
                row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + li
                for lj in 1:p.n_rot
                    if p.kDDmat[li, lj] > 0
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

            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
            val = p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            ### pumping for U level
            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
            col = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + (p.JL-p.K0+1) #
            val = p.pump_IR[vi, ri]
            s = put_row_col_val(rowind, colind, value, row, col, val, s)

            col = row
            val = - p.pump_IR[vi, ri] * p.g_L/p.g_U
            s = put_row_col_val(rowind, colind, value, row, col, val, s)
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
    # first layer, Neumann BC, rotational levels
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
    # Robin BC for rotational levels:
    vbar = p.v_avg/2/sqrt(2)/p.norm_time
    for vi in 1:p.num_freq
        for j in 1:p.n_rot
            row = p.layer_unknown*Ri + (vi-1)*p.n_rot + j
            s = put_row_col_val(rowind, colind, value, row, row, -p.D/p.Δr - vbar/4, s)
            s = put_row_col_val(rowind, colind, value, row, row-p.layer_unknown, p.D/p.Δr - vbar/4, s)
            for lprime in 1:p.n_rot # rot level with the same velocity
                col = p.layer_unknown*Ri + (vi-1)*p.n_rot + lprime
                s = put_row_col_val(rowind, colind, value, row, col, vbar/4 * p.rotpopfr[j], s)
                s = put_row_col_val(rowind, colind, value, row, col-p.layer_unknown, vbar/4 * p.rotpopfr[j], s)
            end
            # flux back from V Sigma into level j and vel vi
            col = p.layer_unknown*(Ri+1)
            s = put_row_col_val(rowind, colind, value, row, col, vbar/4 * p.rotpopfr[j]*p.gauss_dist[vi], s)
            s = put_row_col_val(rowind, colind, value, row, col-p.layer_unknown, vbar/4 * p.rotpopfr[j]*p.gauss_dist[vi], s)
        end
    end

    # for vib Sigma level
    row = p.layer_unknown*(Ri+1)
    s = put_row_col_val(rowind, colind, value, row, row, -p.D/p.Δr - vbar/4, s)
    s = put_row_col_val(rowind, colind, value, row, row-p.layer_unknown, p.D/p.Δr - vbar/4, s)
    # flux back from all rot levels:
    for j in 1:p.layer_unknown
        col = p.layer_unknown*Ri + j
        s = put_row_col_val(rowind, colind, value, row, col, vbar/4 * p.f_6_0, s)
        s = put_row_col_val(rowind, colind, value, row, col-p.layer_unknown, vbar/4 * p.f_6_0, s)
    end

    ########################## rot-vib transition terms for rot levels ##########################

    for ri in 1:p.num_layers
        for vi in 1:p.num_freq
            for li in 1:p.n_rot # V0
                row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + li
                col = ri*p.layer_unknown
                if li > p.n_rot÷2 # in V3
                    s = put_row_col_val(rowind, colind, value, row, row, -A, s)
                    s = put_row_col_val(rowind, colind, value, row, col, A*p.f_3A[ri]/p.f_6A[ri]*p.cj[li]*p.gauss_dist[vi], s)
                else # in V0
                    s = put_row_col_val(rowind, colind, value, row, row, -B*p.f_6A[ri]/p.f_GA[ri], s)
                    s = put_row_col_val(rowind, colind, value, row, row, B*p.cj[li]*p.gauss_dist[vi], s)
                end
            end
        end
    end

    ###################### vib V sigma ##################################
    for ri in 1:p.num_layers
        # V Sigma
        row = ri*p.layer_unknown
        s = put_row_col_val(rowind, colind, value, row, row, 1., s)
        # for vi in 1:p.num_freq
        #     for l in 1:p.n_rot÷2
        #         col = p.layer_unknown*(ri-1) + (vi-1)*p.n_rot + p.n_rot÷2 + l
        #         s = put_row_col_val(rowind, colind, value, row, col, -p.k36A[ri]/p.k63A[ri], s)
        #     end
        # end

        for vi in 1:p.num_freq
            for l in 1:p.n_rot
                col = p.layer_unknown*(ri-1) + (vi-1)*p.n_rot + l
                s = put_row_col_val(rowind, colind, value, row, col, 1.0, s)
            end
        end
    end

    # row = (p.num_layers+1) * p.layer_unknown
    # for j in 1:p.num_layers*p.layer_unknown
    #     ri = (j-1)÷p.layer_unknown + 1
    #     s = put_row_col_val(rowind, colind, value, row, j, p.r_int[ri], s)
    # end

end

function totN0r(p, sol, ri)
    tot = 0.
    for vi in 1:p.num_freq
        rowoffset = p.layer_unknown*(ri-1) + (vi-1)*p.n_rot
        tot += sum(sol[rowoffset+1:rowoffset+p.n_rot÷2])
    end
    return tot + p.ntotal*p.f_G_0
end

function totN3r(p, sol, ri)
    tot = 0.
    for vi in 1:p.num_freq
        rowoffset = p.layer_unknown*(ri-1) + (vi-1)*p.n_rot + p.n_rot÷2
        tot += sum(sol[rowoffset+1:rowoffset+p.n_rot÷2])
    end
    return tot + p.ntotal*p.f_3_0
end
