function pop_inv_dir_layer(para::Param, sol::AbstractVector, layer::Integer)
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib;
    index_offset = layer_unknown * (layer-1);
    pop_inv = 0;
    for vi in 1:para.num_freq
        index_u = index_offset + (vi-1) * n_rot + 12;
        pop_inv = pop_inv + (sol[index_u]/g_U - sol[index_u - 1]/g_L) * para.fp_lasing[vi];
    end
    return pop_inv
end

function get_Nu_vs_freq_layer(para::parameter, sol::Array, layer::Int64)
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib;
    index_offset = layer_unknown * (layer-1);
    Nu = zeros(para.num_freq)
    for i in 1:para.num_freq
        index_u = index_offset + (i-1) * n_rot + 12;
        Nu[i] = sol[index_u]
    end
    return Nu
end

function get_total_Nu_dist_layer(para, sol, layer)
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib
    index_p = layer_unknown * layer - n_vib + 2
    Nu_sol = get_Nu_vs_freq_layer(para, sol, layer)
    Nu = Nu_sol + C5U * (para.ntotal * f_3/2 + sol[index_p]) .* para.gauss_dist
    return Nu
end

function get_Nl_vs_freq_layer(para::parameter, sol::Array, layer::Int64)
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib;
    index_offset = layer_unknown * (layer-1);
    Nl = zeros(para.num_freq)
    for i in 1:para.num_freq
        index_l = index_offset + (i-1) * n_rot + 2;
        Nl[i] = sol[index_l]
    end
    return Nl
end

function get_total_Nl_dist_layer(para, sol, layer)
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib
    index_p = layer_unknown * layer - n_vib + 1
    Nl_sol = get_Nl_vs_freq_layer(para, sol, layer)
    Nl = Nl_sol + (para.ntotal * f_G/2 + sol[index_p]) .* para.gauss_dist
    return Nl
end

function get_rot_vs_freq_layer(para::parameter, sol::Array, layer::Int64, level::Int64)
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib;
    index_offset = layer_unknown * (layer-1);
    N = zeros(para.num_freq)
    for i in 1:para.num_freq
        index = index_offset + (i-1) * n_rot + level;
        N[i] = sol[index]
    end
    return N
end

function get_inv_U(para::parameter, sol::Array, layer::Int64)
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib
    Nu_sol = get_Nu_vs_freq_layer(para, sol, layer)
    Nu_1_sol = get_rot_vs_freq_layer(para, sol, layer, 11)
    index_p = layer_unknown * layer - n_vib + 2
    Nu = get_total_Nu_dist_layer(para, sol, layer)
    Nu_1 = Nu_1_sol + C4U * (para.ntotal * f_3/2 + sol[index_p]) .* para.gauss_dist
    return (Nu - Nu_1*g_U/g_L)
end

function get_inv_L(para::parameter, sol::Array, layer::Int64)
    layer_unknown::Int64 = n_rot*para.num_freq + n_vib
    Nl_sol = get_Nl_vs_freq_layer(para, sol, layer)
    Nl_1_sol = get_rot_vs_freq_layer(para, sol, layer, 3)
    index_p = layer_unknown * layer - n_vib + 1
    Nl = get_total_Nl_dist_layer(para, sol, layer)
    Nl_1 = Nl_1_sol + C5L * (para.ntotal * f_G/2 + sol[index_p]) .* para.gauss_dist
    return (Nl_1 - Nl*g_U/g_L)
end
