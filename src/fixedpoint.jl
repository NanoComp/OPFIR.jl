function fixedpoint(sol_0, p, matrix_0, lu_mat0)
    max_ele = p.num_freq * p.num_layers * (p.n_rot*(p.n_rot+2) + p.n_vib*(p.n_rot+p.n_vib+2))*2

    rowind = ones(Int64, max_ele)
    colind = ones(Int64, max_ele)
    value = zeros(max_ele)
    rhs = zeros(p.num_layers*p.layer_unknown + p.layer_unknown)

    update_alpha_from_N!(p, sol_0)
    update_Param_from_alpha!(p, sol_0)
    println("absorption: ", p.alpha_r[1])

    compute_rhs(rhs, p, sol_0)
    compute_row_col_val(rowind, colind, value, p, sol_0)
    println("finished computing row, col, value!")
    matrix = sparse(rowind, colind, value)
    mat_rhs_modify(matrix, rhs, p)

    println("start to compute the solution!")
    sol_1 = matrix \ rhs
    # println(matrix[p.layer_unknown-p.n_vib+1, :], ", rhs:", rhs[p.layer_unknown-p.n_vib+1])

    update_alpha_from_N!(p, sol_1)

    if p.optcavity && p.WiU == 0. && p.WiL == 0.
        p.L = 0.5/p.alpha_r[1]*100
        println("L = ", p.L, "cm")
    else
        println("L = ", p.L, "cm")
    end

    println("norm of sol diff = ", norm(sol_1 - sol_0) / norm(sol_1))
    flush(STDOUT)

    return sol_1
end

function mat_rhs_modify(matrix, rhs, p)
    for ri in p.num_layers
        row = ri*p.layer_unknown
        matrix[row,:] = 0.0
        rhs[row] = 0.0
        for col in (ri-1)*p.layer_unknown + 1 : ri*p.layer_unknown
            matrix[row, col] = 1.0
        end
    end
end
