function fixedpoint(sol_0, p)

    temp_flag = p.solstart_flag
    # if p.solstart_flag==0
    # else
    #     tot_pump = 0.0
    #     # p.T_vA = Tv(p, sol_0[end-p.n_vib+1]+p.ntotal*p.f_G_0/2,
    #     #              sol_0[end-p.n_vib+2]+p.ntotal*p.f_3_0/2)
    #     # p.T_vA = 300
    #     Q_v = Qv_0(p, p.T_vA)
    #     p.f_G[:] = fraction_Vg(p, p.T_vA)
    #     p.f_3[:] = fraction_V3(p, p.T_vA)
    #     p.f_6[:] = 1 - p.f_G[:] - p.f_3[:]
    #     # for ri in 1:p.num_layers
    #     #     # pump_ri = pump_total(p, sol_0, ri)*p.f_6[ri]/(p.f_3[ri]+p.f_6[ri])
    #     #     # Q_v = p.Q/(1-pump_total(p, sol_0, ri)*Q_v/p.kwall[ri]/p.ntotal)
    #     #     # if kwall very small, Q_v -> negative, not possible
    #     #     # T_v = solve_Tv(Q_v, p)
    #     #     p.netrate_36A[ri] = pump_total(p, sol_0, ri) * p.f_6[ri]*Q_v/(Q_v-1)
    #     # end
    #     p.k36 = 100.0
    #     p.k63 = p.k36 * p.f_3[1]/p.f_6[1]
    #     println("Tv = ", p.T_vA, " k36 = ", p.k36, " k63 = ", p.k63)
    #     p.solstart_flag=0
    # end
    # Q_v = Qv_0(p, p.T_vA)
    # if p.T_vA != 300
        p.f_G[:] = fraction_Vg(p, p.T_vA)
        p.f_3[:] = fraction_V3(p, p.T_vA)
        p.f_6[:] = 1 - p.f_G[:] - p.f_3[:]
        p.k36 = 100.0
        p.k63 = p.k36 * p.f_3[1]/p.f_6[1]
    # end

    println("k36 = ", p.k36, " k63 = ", p.k63)

    max_ele = p.num_freq * p.num_layers * (p.n_rot*(p.n_rot+2) + p.n_vib*(p.n_rot+p.n_vib+2))

    rowind = ones(Int64, max_ele)
    colind = ones(Int64, max_ele)
    value = zeros(max_ele)
    rhs = zeros(p.num_layers*p.layer_unknown)

    compute_rhs(rhs, p, sol_0)
    compute_row_col_val(rowind, colind, value, p, sol_0)

    matrix = sparse(rowind, colind, value)

    for j = 1:size(matrix,2)
        matrix[1, j] = 0
        matrix[end, j] = 0
    end

    for j = 1:p.num_layers
        index1 = p.layer_unknown * (j-1) + 1
        index2 = p.layer_unknown * j
        matrix[1, index1:index2] = p.r_int[j]

        index1 = p.layer_unknown * j - p.n_vib + p.n_vib÷2 + 1
        index2 = p.layer_unknown * j
        matrix[end, index1:index2] = p.r_int[j]
    end

    rhs[1] = 0
    rhs[end] = 0
    if p.lin_solver=="Sherman_Morrison"
        sol_1 = - SM_solve(matrix, 2, [1, size(matrix,1)], rhs)
    elseif p.lin_solver=="Default"
        sol_1 = - matrix \ rhs
    end

    println("norm of sol diff = ", norm(sol_1 - sol_0) / norm(sol_1))
    p.solstart_flag = temp_flag
    return sol_1
end

function update_f(p, sol)

end
