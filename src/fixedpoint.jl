function fixedpoint(sol_0, p, matrix_0, lu_mat0)
    # max_ele = p.num_freq * p.num_layers * (p.n_rot*(p.n_rot+2) + p.n_vib*(p.n_rot^2+p.n_vib+2))

    # rowind = ones(Int64, max_ele)
    # colind = ones(Int64, max_ele)
    # value = zeros(max_ele)
    rowind = Array{Float64}(0)
    colind = Array{Float64}(0)
    value = Array{Float64}(0)
    rhs = zeros((p.num_layers+1)*p.layer_unknown)

    update_alpha_from_N!(p, sol_0)
    update_Param_from_alpha!(p, sol_0)

    # if p.model_flag == 2
    #     status = 1
    #     status = updateTv(p, sol_0)
    #     if status == 0
    #         return 0
    #     end
    #     updateks(p)
    # end
    compute_rhs(rhs, p, sol_0)
    compute_row_col_val(rowind, colind, value, p, sol_0)

    matrix = sparse(rowind, colind, value)
    if p.model_flag == 1
        mat_rhs_modify(matrix, rhs, p)
    end

    if p.solstart_flag == 1
        matrix_B0 = matrix - matrix_0
        rhs = rhs - matrix_B0 * sol_0
        sol_1 = lu_mat0 \ rhs
    else
        sol_1 = matrix \ rhs
    end
    # println(matrix[p.layer_unknown-p.n_vib+1, :], ", rhs:", rhs[p.layer_unknown-p.n_vib+1])

    update_alpha_from_N!(p, sol_1)

    if p.optcavity && p.WiU == 0. && p.WiL == 0.
        p.L = 0.4/p.alpha_r[1]*100
        println("L = ", p.L, "cm")
    else
        println("L = ", p.L, "cm")
    end
    println(size(matrix))
    println("length of sol1: ", length(sol_1))
    println("norm of sol diff = ", norm(sol_1 - sol_0) / norm(sol_1))
    flush(STDOUT)

    return sol_1
end

function mat_modify(matrix, p)
    for ri in 1:p.num_layers
        # V26A and V26E:
        rowA = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 6
        for row in vcat(rowA, rowA + p.n_vib÷2)
            for k in 0:5
                matrix[row, row-k] = 1.
            end
        end
    end
end

function mat_rhs_modify(matrix, rhs, p)
    mat_modify(matrix, p)
    # rhs[1] = 0
    # rhs[end] = 0
end
