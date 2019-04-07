function fixedpoint(sol_0, p, matrix_0, lu_mat0)
    max_ele = p.num_freq * p.num_layers * (p.n_rot*(p.n_rot+2) + p.n_vib*(p.n_rot+p.n_vib+2))*2

    rowind = ones(Int64, max_ele)
    colind = ones(Int64, max_ele)
    value = zeros(max_ele)
    rhs = zeros(p.num_layers*p.layer_unknown + p.layer_unknown)

    update_alpha_from_N!(p, sol_0)
    update_Param_from_alpha!(p, sol_0)
    # println("absorption: ", p.alpha_r[1])

    compute_rhs(rhs, p, sol_0)

    compute_row_col_val(rowind, colind, value, p, sol_0)

    matrix = sparse(rowind, colind, value)
    mat_rhs_modify(matrix, rhs, p)
    # mat_rhs_modify(rowind, colind, value, rhs, p)
    # println("start to compute the solution!")

    # sol_1 = matrix \ rhs
    sol_1 = solveMUMPS(matrix, rhs)
    # println(matrix[p.layer_unknown-p.n_vib+1, :], ", rhs:", rhs[p.layer_unknown-p.n_vib+1])

    update_alpha_from_N!(p, sol_1)

    if p.optcavity && p.WiU == 0. && p.WiL == 0.
        p.L = 0.5/p.alpha_r[1]*100
        println("L = ", p.L, "cm")
    else
        # println("L = ", p.L, "cm")
    end

    println("norm of sol diff = ", norm(sol_1 - sol_0) / norm(sol_1))
    flush(STDOUT)

    return sol_1
end

function mat_rhs_modify(matrix, rhs, p)
    for ri in 1:p.num_layers+1
        row = (ri-1)*p.layer_unknown + 1
        # for vi in 1:p.num_freq
        #     row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 1
        matrix[row,:] = 0.0
        rhs[row] = 0.0
        # for col in row : row+p.layer_unknown-1
        matrix[row, row : row+p.layer_unknown-1] = 0.1
        # end
        # end
    end
end

# function mat_rhs_modify(rowind, colind, value, rhs, p)
#     for ri in 1:p.num_layers
#         for vi in 1:p.num_freq
#             row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 1
#             rhs[row] = 0.0
#             for r in 1:length(rowind)
#                 if rowind[r] == row
#                     value[r] = 0.0
#                 end
#             end
#             for col in row+1:row+p.n_rot
#                 push!(rowind, row)
#                 push!(colind, col)
#                 push!(value, 1.0)
#             end
#         end
#     end
# end
