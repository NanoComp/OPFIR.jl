using MUMPSjInv

function fixedpoint(sol_0, p, matrix_0, lu_mat0; mumps_solver=0)
    # max_ele = p.num_freq * p.num_layers * (p.n_rot*(p.n_rot+2) + p.n_vib*(p.n_rot^2+p.n_vib+2))

    # rowind = ones(Int64, max_ele)
    # colind = ones(Int64, max_ele)
    # value = zeros(max_ele)

    rowind = Array{Float64}(0)
    colind = Array{Float64}(0)
    value = Array{Float64}(0)
    rhs = zeros(p.num_layers*p.layer_unknown + p.n_vib)

    update_alpha_from_N!(p, sol_0)
    update_Param_from_alpha!(p, sol_0)

    # tic()
    compute_rhs(rhs, p, sol_0)
    compute_row_col_val(rowind, colind, value, p, sol_0)
    # toc()
    matrix = sparse(rowind, colind, value)
    mat_rhs_modify(matrix, rhs, p)
    tic()
    if p.solstart_flag == 1
        matrix_B0 = matrix - matrix_0
        rhs = rhs - matrix_B0 * sol_0
        sol_1 = lu_mat0 \ rhs
    else
        if mumps_solver == 0
            sol_1 = matrix \ rhs
        elseif mumps_solver == 1
            sol_1 = solveMUMPS(matrix, rhs)
        elseif mumps_solver == 2
            sol_1 = MUMPS.solve(SparseMatrixCSC{Float64,Int32}(matrix), rhs)
        end
    end
    # println(matrix[p.layer_unknown-p.n_vib+1, :], ", rhs:", rhs[p.layer_unknown-p.n_vib+1])

    update_alpha_from_N!(p, sol_1)
    toc()
    if p.optcavity && p.WiU == 0. && p.WiL == 0.
        p.L = 0.4/p.alpha_r[1]*100
        println("L = ", p.L, "cm")
    else
        println("L = ", p.L, "cm")
    end

    println("norm of sol diff = ", norm(sol_1 - sol_0) / norm(sol_1))
    flush(STDOUT)

    return sol_1
end

function mat_modify(matrix, p)
    for ri in 1:p.num_layers
        row = ri*p.layer_unknown
        # matrix[row, :] = 0
        for k in vcat(0:p.n_vib-1)
            matrix[row, row-k] = 0.1
        end
    end
end

function mat_rhs_modify(matrix, rhs, p)
    mat_modify(matrix, p)
    # rhs[1] = 0
    # rhs[end] = 0
end
