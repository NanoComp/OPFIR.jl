function fixedpoint(sol_0, p)
    # if initial_flag==0
    # else
    #     tot_pump = 0.0
    #     temp_36 = 0.0
    #     temp_63 = 0.0
    #     for ri in 1:p.num_layers
    #         tot_pump += pump_total(p, sol_0, ri)*ri
    #         ind3 = p.layer_unknown*ri - p.n_vib + 2
    #         temp_36 += (sol_0[ind3] - p.f_3[ri]/p.f_6[ri] * sol_0[ind3+1])*ri
    #         temp_63 += (p.f_6[ri]/p.f_3[ri] * sol_0[ind3] - sol_0[ind3+1])*ri
    #     end
    #     p.k63 = tot_pump/(temp_36+1e-6)
    #     p.k36 = p.f_6_0/p.f_3_0 * p.k63
    # end
    # p.k36 = tot_pump/(temp_63+1e-6)

    # println("k36 = ", p.k36, " k63 = ", p.k63)

    max_ele = p.num_freq * p.num_layers * (p.n_rot*(p.n_rot+2) + p.n_vib*(p.n_rot+p.n_vib+2))

    rowind = ones(Int64, max_ele)
    colind = ones(Int64, max_ele)
    value = zeros(max_ele)
    rhs = zeros(p.num_layers*p.layer_unknown)

    if p.model_flag==2
        update_f(p, sol_0)
    end
    compute_rhs(rhs, p, sol_0)
    compute_row_col_val(rowind, colind, value, p, sol_0)
    # if p.model_flag==2
    # println(p.ntotal, sol_0[p.layer_unknown-p.n_vib+1:p.layer_unknown])
    # end

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

    return sol_1
end

function update_f(p, sol)

end
