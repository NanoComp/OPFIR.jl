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

# population inversion of direct lasing in each radial layer
function pop_inv_dir_layer(p, sol, layer)
    index_offset = ind_offset(p, layer)
    pop_inv = 0
    for vi in 1:p.num_freq
        index_u = index_offset + (vi-1) * p.n_rot + 12
        pop_inv += (sol[index_u]/p.g_U - sol[index_u - 1]/p.g_L) * p.fp_lasing[vi]
    end
    return pop_inv
end

# nonthermal/thermal/total Nu distribution (over freq/velocity) in each radial layer
function Nu_NT_dist_layer(p, sol, layer)
    index_offset = ind_offset(p, layer)
    Nu = zeros(p.num_freq)
    for i in 1:p.num_freq
        index_u = index_offset + (i-1) * p.n_rot + 12
        Nu[i] = sol[index_u]
    end
    return Nu
end
function Nu_T_dist_layer(p, sol, layer)
    index_p = ind_p_U(p, layer)
    return p.C5U * (p.ntotal * p.f_3_0/2 + sol[index_p]) .* p.gauss_dist
end

function Nu_total_dist_layer(p, sol, layer)
    return Nu_NT_dist_layer(p, sol, layer) + Nu_T_dist_layer(p, sol, layer)
end

# nonthermal/thermal/total Nu-1 distribution (over freq/velocity) in each radial layer
function Nu_1_NT_dist_layer(p, sol, layer)
    index_offset = ind_offset(p, layer)
    N = zeros(p.num_freq)
    for i in 1:p.num_freq
        index_u = index_offset + (i-1) * p.n_rot + 11
        N[i] = sol[index_u]
    end
    return N
end
function Nu_1_T_dist_layer(p, sol, layer)
    index_p = ind_p_U(p, layer)
    return p.C4U * (p.ntotal * p.f_3_0/2 + sol[index_p]) .* p.gauss_dist
end

function Nu_1_total_dist_layer(p, sol, layer)
    return Nu_1_NT_dist_layer(p, sol, layer) + Nu_1_T_dist_layer(p, sol, layer)
end

# nonthermal/thermal/total Nl distribution (over freq/velocity) in each radial layer
function Nl_NT_dist_layer(p, sol, layer)
    index_offset = ind_offset(p, layer)
    Nl = zeros(p.num_freq)
    for i in 1:p.num_freq
        index_l = index_offset + (i-1) * p.n_rot + 2
        Nl[i] = sol[index_l]
    end
    return Nl
end
function Nl_T_dist_layer(p, sol, layer)
    index_p = ind_p_L(p, layer)
    return p.C4L * (p.ntotal * p.f_G_0/2 + sol[index_p]) .* p.gauss_dist
end

function Nl_total_dist_layer(p, sol, layer)
    return Nl_NT_dist_layer(p, sol, layer) + Nl_T_dist_layer(p, sol, layer)
end

# nonthermal/thermal/total Nl+1 distribution (over freq/velocity) in each radial layer
function Nl_1_NT_dist_layer(p, sol, layer)
    index_offset = ind_offset(p, layer)
    N = zeros(p.num_freq)
    for i in 1:p.num_freq
        index_l = index_offset + (i-1) * p.n_rot + 3
        N[i] = sol[index_l]
    end
    return N
end
function Nl_1_T_dist_layer(p, sol, layer)
    index_p = ind_p_L(p, layer)
    return p.C5L * (p.ntotal * p.f_G_0/2 + sol[index_p]) .* p.gauss_dist
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
    return (Nu - Nu_1*p.g_U/p.g_L)
end

# pop inv distribution over velocity between L+1 and L
function inv_L_dist_layer(p, sol, layer)
    Nl = Nl_total_dist_layer(p, sol, layer)
    Nl_1 = Nl_1_total_dist_layer(p, sol, layer)
    return (Nl_1 - Nl*p.g_U/p.g_L)
end

##
function gain_dir_layer(p, sol, layer)
    inv_U = inv_U_dist_layer(p, sol, layer)
    inv = 0.0
    for i in 1:length(inv_U)
#        inv += inv_U[i]
        g(ν) = 1/π * p.Δ_fP ./ ((ν - p.f_dist_dir_lasing[i]).^2 + p.Δ_fP^2)
#        fraction = quadgk(g, p.f_dir_lasing-3e6, p.f_dir_lasing+3e6)[1]
        fraction = 1
        inv += inv_U[i] * fraction
    end
##    sum(inv_U .* p.fp_lasing)
    return inv
end

function gain_ref_layer(p, sol, layer)
    inv_L = inv_L_dist_layer(p, sol, layer)
    inv = 0.0
    for i in 1:length(inv_L)
        # inv += inv_L[i]
        g(ν) = 1/π * p.Δ_fP ./ ((ν - p.f_dist_ref_lasing[i]).^2 + p.Δ_fP^2)
 #       fraction = quadgk(g, p.f_ref_lasing-3e6, p.f_ref_lasing+3e6)[1]
         fraction = 1.0
         inv += inv_L[i] * fraction
    end
    return inv
end


function gain_profile_layer(p, sol, layer; ω=p.f_dist_dir_lasing)
    dist = Array(Float64, length(ω))
    inv_U = inv_U_dist_layer(p, sol, layer)
    for i in 1:length(ω)
        ωi = ω[i]
        ω_dist = f_NT_ampl(ω, p.Δ_f_NT, ωi)
        ω_dist = ω_dist/sum(ω_dist)
        dist[i] = sum(inv_U .* ω_dist)
    end
    return dist
end

### gain calculation with integration over the field: ###
function mode(mode_num)
    if mode_num == 1 || mode_num == 3 || mode_num == 7
        m = 0
    elseif mode_num == 2 || mode_num == 5 || mode_num == 8
        m = 1
    elseif mode_num == 4 || mode_num == 6
        m = 2
    end
    return m
end

function gain_dir(p, sol; LasLevel="U")
    m = mode(p.mode_num)
    radius_m = p.radius/100
    k_rol_library = p.p_library/radius_m
    k_rol = k_rol_library[p.mode_num]

    ϕ_list = linspace(0, 2*pi, 150)
    denom = 0.0
    numerator = 0.0
    pop_inversion = zeros(size(p.r_int))
    gain_LU = 0.0

    for ri in 1:p.num_layers
        r = p.r_int[ri]
        if LasLevel=="U"
            pop_inversion[ri] = gain_dir_layer(p, sol, ri)
        elseif LasLevel=="L"
            pop_inversion[ri] = gain_ref_layer(p, sol, ri)
        else
            println("Only support L or U lasing!")
        end

        for ϕ in ϕ_list
            if p.mode_num >= 7 ## TM modes
                Ez = besselj(m, k_rol*r) * cos(m*ϕ)
                Er = derv_bessel(m, k_rol*r) * cos(m*ϕ)
                Eϕ = m/(k_rol*r) * besselj(m, k_rol*r) * (-sin(m*ϕ))
            else ## TE modes
                Ez = 0
                Er = m/(k_rol*r) * besselj(m, k_rol*r) * (-sin(m*ϕ))
                Eϕ = derv_bessel(m, k_rol*r) * cos(m*ϕ)
            end
            E_sq = Ez*conj(Ez) + Er*conj(Er) + Eϕ*conj(Eϕ)
            denom += E_sq * r # cylindrical coordinate
            numerator += E_sq * r * pop_inversion[ri]
        end
    end
    if LasLevel=="U"
        λ_las = p.c/p.f_dir_lasing
    elseif LasLevel=="L"
        λ_las = p.c/p.f_ref_lasing
    end
    gain_LU = numerator/denom *
           λ_las^2/(8*pi*p.n0^2*p.t_spont)/p.Δν_THz*0.01 # in cm^-1
    return gain_LU #* (p.L_eff/p.L)
end

function derv_bessel(ν,x)
    return 0.5.*(besselj(ν-1,x)-besselj(ν+1,x))
end
