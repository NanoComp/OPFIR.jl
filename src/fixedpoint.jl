function fixedpoint(sol_0, p, matrix_0, lu_mat0)
    max_ele = p.num_freq * p.num_layers * (p.n_rot*(p.n_rot+2) + p.n_vib*(p.n_rot+p.n_vib+2))

    rowind = ones(Int64, max_ele)
    colind = ones(Int64, max_ele)
    value = zeros(max_ele)
    rhs = zeros(p.num_layers*p.layer_unknown)

    update_alpha_from_N!(p, sol_0)
    # println(p.alpha_r)
    update_Param_from_alpha!(p, sol_0)

    if p.model_flag == 2
        updateTv(p, sol_0)
        updateks(p)
    end
    # println(p.alpha_r)
    # println(p_0.alpha_r)
    compute_rhs(rhs, p, sol_0)
    compute_row_col_val(rowind, colind, value, p, sol_0)

    matrix = sparse(rowind, colind, value)
    mat_rhs_modify(matrix, rhs, p)

    if p.solstart_flag == 1
        # rowind_0 = ones(Int64, max_ele)
        # colind_0 = ones(Int64, max_ele)
        # value_0 = zeros(max_ele)
        #
        # compute_row_col_val(rowind_0, colind_0, value_0, p_0, sol_in)
        # matrix_0 = sparse(rowind_0, colind_0, value_0)
        # mat_modify(matrix_0, p_0)
        matrix_B0 = matrix - matrix_0
        # println(norm(matrix_B0, 1), " ", norm(matrix_0, 1), " ", norm(matrix,1))
        # sol_1 = - matrix \ rhs
        rhs = rhs + matrix_B0 * sol_0
        # mat_rhs_modify(matrix_0, rhs, p_0)
        sol_1 = lu_mat0 \ (-rhs)

    else
        sol_1 = - matrix \ rhs
    end
    # sol_1 = - matrix \ rhs
    # if p.lin_solver=="Sherman_Morrison"
    #     sol_1 = - SM_solve(matrix, 2, [1, size(matrix,1)], rhs)
    # elseif p.lin_solver=="Default"
    #     sol_1 = - matrix \ rhs
    # end
    update_alpha_from_N!(p, sol_1)

    # println(p.averagePF[1], ", ", p.averagePF[end])
    optimizecavity = 0
    if optimizecavity == 1 && p.WiU == 0. && p.WiL == 0.
        p.L = 1/p.alpha_r[1]*100
        println("L = ", p.L, "cm")
    else
        println("L = ", p.L, "cm")
    end
    #p.WiU!=0. && throw(ArgumentError("wi not equal to zero, L should be fixed!"))

    println("norm of sol diff = ", norm(sol_1 - sol_0) / norm(sol_1))
    flush(STDOUT)

    return sol_1
end

function mat_modify(matrix, p)
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
end

function mat_rhs_modify(matrix, rhs, p)
    mat_modify(matrix, p)
    rhs[1] = 0
    rhs[end] = 0
end
