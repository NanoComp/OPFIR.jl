using MUMPSjInv

function fixedpoint(sol_0, p, matrix_0, lu_mat0)
    rowind = Array{Float64}(undef, 0)
    colind = Array{Float64}(undef, 0)
    value = Array{Float64}(undef, 0)
    rhs = zeros(p.num_layers*p.layer_unknown + p.n_vib)

    update_alpha_from_N!(p, sol_0)
    update_Param_from_alpha!(p, sol_0)

    compute_rhs(rhs, p, sol_0)
    compute_row_col_val(rowind, colind, value, p, sol_0)

    matrix = sparse(rowind, colind, value)
    mat_rhs_modify(matrix, rhs, p)

    if p.solstart_flag == 1
        matrix_B0 = matrix - matrix_0
        rhs = rhs - matrix_B0 * sol_0
        sol_1 = lu_mat0 \ rhs
    else
        if p.mumps_solver == 0
            sol_1 = matrix \ rhs
        elseif p.mumps_solver == 1
            sol_1 = solveMUMPS(matrix, rhs)
        elseif p.mumps_solver == 2
            sol_1 = MUMPS.solve(SparseMatrixCSC{Float64,Int32}(matrix), rhs)
        end
    end

    update_alpha_from_N!(p, sol_1)

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
