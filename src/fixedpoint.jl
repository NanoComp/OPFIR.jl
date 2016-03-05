function fixedpoint(sol_0, p)
    temp_flag = p.solstart_flag

#    println("k36 = ", p.k36A[1:2], " k63 = ", p.k63A[1:2])

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
    flush(STDOUT)
    p.solstart_flag = temp_flag
    return sol_1
end
