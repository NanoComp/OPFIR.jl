function OPFIR_fixedpoint(sol_0, para)
    layer_unknown = n_rot*para.num_freq + n_vib
    max_ele = para.num_freq * para.num_layers * (n_rot*(n_rot+2) + n_vib*(n_rot+n_vib+2))

    rowind = ones(Int64,max_ele)
    colind = ones(Int64,max_ele)
    value = zeros(max_ele)
    rhs = zeros(para.num_layers*layer_unknown)

    OPFIR_compute_rhs(rhs, para)
    OPFIR_compute_row_col_val(rowind, colind, value, para, sol_0)

    matrix = sparse(rowind, colind, value)

    ############################## deal with the singularity ####################
    for j = 1:size(matrix,1)
        matrix[1,j] = 0;
        matrix[end,j] = 0;
    end

    for j = 1:para.num_layers
        index1::Int64 = layer_unknown * (j-1) + 1
        index2::Int64 = layer_unknown * j
        matrix[1, index1:index2] = para.r_int[j]

        # ind_0E::Int64 = layer_unknown * j - n_vib + 7;
        # ind_3E::Int64 = layer_unknown * j - n_vib + 8;
        ## close the two E type levels:
        # matrix[index2-1, index2-1] = 1.0;
        # matrix[index2-1, index2] = 1.0;
        index1 = layer_unknown * j - n_vib + 7
        index2 = layer_unknown * j - n_vib + 12
        matrix[end, index1:index2] = para.r_int[j]
    end

    rhs[1] = 0
    rhs[end] = 0

    sol_1 = - matrix \ rhs
    # sol_1 = - SM_solve(matrix, 2, [1, size(matrix,1)], rhs)
    println("norm of sol diff = ", norm(sol_1 - sol_0) / norm(sol_1))
    # println("peak of Nu(v) = ", maximum(sol_0[(12:n_rot:n_rot*para.num_freq)]))
    return sol_1
end
