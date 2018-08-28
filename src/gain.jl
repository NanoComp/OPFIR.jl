# index offset for radial layer n
function ind_offset(p, n)
    return p.layer_unknown * (n -1)
end
# index of the thermal pool U and L are in, i.e., V3 and V0
function ind_p_U(p, layer)
    return p.layer_unknown * layer - p.n_vib + 2
end

function ind_p_L(p, layer)
    return p.layer_unknown * layer - p.n_vib + 1
end

# nonthermal/thermal/total Nu distribution (over freq/velocity) in each radial layer
function Nu_NT_dist_layer(p, sol, layer)
    index_offset = ind_offset(p, layer)
    Nu = zeros(p.num_freq)
    for i in 1:p.num_freq
        index_u = index_offset + (i-1) * p.n_rot + (p.n_rot÷2 + p.JU-p.K0+1)
        Nu[i] = sol[index_u]
    end
    return Nu
end
function Nu_T_dist_layer(p, sol, layer)
    index_p = ind_p_U(p, layer)
    return p.CU * (p.ntotal * p.f_3_0 + sol[index_p]) .* p.gauss_dist
    return p.CU * (p.ntotal * p.f_3_0/2 + sol[index_p]) .* p.gauss_dist
end

function Nu_total_dist_layer(p, sol, layer)
    return Nu_NT_dist_layer(p, sol, layer) + Nu_T_dist_layer(p, sol, layer)
end

# nonthermal/thermal/total Nu-1 distribution (over freq/velocity) in each radial layer
function Nu_1_NT_dist_layer(p, sol, layer)
    index_offset = ind_offset(p, layer)
    N = zeros(p.num_freq)
    for i in 1:p.num_freq
        index_u = index_offset + (i-1) * p.n_rot + (p.n_rot÷2 + p.JU-p.K0)
        N[i] = sol[index_u]
    end
    return N
end
function Nu_1_T_dist_layer(p, sol, layer)
    index_p = ind_p_U(p, layer)
    return p.CU1 * (p.ntotal * p.f_3_0 + sol[index_p]) .* p.gauss_dist
    return p.CU1 * (p.ntotal * p.f_3_0/2 + sol[index_p]) .* p.gauss_dist
end

function Nu_1_total_dist_layer(p, sol, layer)
    return Nu_1_NT_dist_layer(p, sol, layer) + Nu_1_T_dist_layer(p, sol, layer)
end

# nonthermal/thermal/total Nl distribution (over freq/velocity) in each radial layer
function Nl_NT_dist_layer(p, sol, layer)
    index_offset = ind_offset(p, layer)
    Nl = zeros(p.num_freq)
    for i in 1:p.num_freq
        index_l = index_offset + (i-1) * p.n_rot + (p.JL-p.K0+1)
        Nl[i] = sol[index_l]
    end
    return Nl
end
function Nl_T_dist_layer(p, sol, layer)
    index_p = ind_p_L(p, layer)
    return p.CL * (p.ntotal * p.f_G_0 + sol[index_p]) .* p.gauss_dist
    return p.CL * (p.ntotal * p.f_G_0/2 + sol[index_p]) .* p.gauss_dist
end

function Nl_total_dist_layer(p, sol, layer)
    return Nl_NT_dist_layer(p, sol, layer) + Nl_T_dist_layer(p, sol, layer)
end

# nonthermal/thermal/total Nl+1 distribution (over freq/velocity) in each radial layer
function Nl_1_NT_dist_layer(p, sol, layer)
    index_offset = ind_offset(p, layer)
    N = zeros(p.num_freq)
    for i in 1:p.num_freq
        index_l = index_offset + (i-1) * p.n_rot + (p.JL-p.K0+2)
        N[i] = sol[index_l]
    end
    return N
end
function Nl_1_T_dist_layer(p, sol, layer)
    index_p = ind_p_L(p, layer)
    return p.CL1 * (p.ntotal * p.f_G_0 + sol[index_p]) .* p.gauss_dist
    return p.CL1 * (p.ntotal * p.f_G_0/2 + sol[index_p]) .* p.gauss_dist
end

function Nl_1_total_dist_layer(p, sol, layer)
    return Nl_1_NT_dist_layer(p, sol, layer) + Nl_1_T_dist_layer(p, sol, layer)
end

# nonthermal distribution of rotational level in each radial layer
function rot_NT_dist_layer(p, sol, layer, level)
    index_offset = ind_offset(p, layer)
    N = zeros(p.num_freq)
    for i in 1:p.num_freq
        index = index_offset + (i-1) * p.n_rot + level
        N[i] = sol[index]
    end
    return N
end

# population inversion distribution over velocity between U and U-1,
function inv_U_dist_layer(p, sol, layer)
    Nu = Nu_total_dist_layer(p, sol, layer)
    Nu_1 = Nu_1_total_dist_layer(p, sol, layer)
    return (Nu - Nu_1*p.g_U/(p.g_U-2))
end

# pop inv distribution over velocity between L+1 and L
function inv_L_dist_layer(p, sol, layer)
    Nl = Nl_total_dist_layer(p, sol, layer)
    Nl_1 = Nl_1_total_dist_layer(p, sol, layer)
    return (Nl_1 - Nl*(p.g_L+2)/p.g_L)
end
