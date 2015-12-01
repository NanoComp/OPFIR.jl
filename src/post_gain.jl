function gain_direct_layer(para::parameter, sol::Array, layer::Int64)
    inv_U = get_inv_U(para, sol, layer) # inversion as a function of velocity
    inv = sum(inv_U .* para.fp_lasing)
    return inv
end

function mode(mode_num)
    if mode_num == 1 || mode_num == 3 || mode_num == 7
        m = 0;
    elseif mode_num == 2 || mode_num == 5 || mode_num == 8
        m = 1;
    elseif mode_num == 4 || mode_num == 6
        m = 2;
    end
    return m
end

function OPFIR_gain(para::parameter, sol::Array)
    m = mode(mode_num)
    radius_m = radius/100
    k_rol_library = p_library/radius_m
    k_rol = k_rol_library[mode_num]

    ϕ_list = linspace(0, 2*pi, 150)
    denom = 0.0
    numerator = 0.0
    pop_inversion = zeros(size(para.r_int))

    for ri::Int64 in 1:para.num_layers
        r = para.r_int[ri]
        pop_inversion[ri] = gain_direct_layer(para, sol, ri)
        for ϕ in ϕ_list
            if mode_num >= 7
                # display('Calculating TM modes....')
                Ez = besselj(m, k_rol*r) * cos(m*ϕ)
                Er = derv_bessel(m, k_rol*r) * cos(m*ϕ)
                Eϕ = m/(k_rol*r) * besselj(m, k_rol*r) * (-sin(m*ϕ))
            else
                # display('Calculating TE modes....');
                Ez = 0
                Er = m/(k_rol*r) * besselj(m, k_rol*r) * (-sin(m*ϕ))
                Eϕ = derv_bessel(m, k_rol*r) * cos(m*ϕ)
            end
            E_sq = Ez*conj(Ez) + Er*conj(Er) + Eϕ*conj(Eϕ)
            denom += E_sq * r
            numerator += E_sq * r * pop_inversion[ri]
        end
    end
    gain = numerator/denom * (c/f_dir_lasing)^2/(8*pi*n0^2*t_spont)/Δν_THz*0.01 # in cm^-1
    return gain
end

function derv_bessel(nu,X)
    return 0.5.*(besselj(nu-1,X)-besselj(nu+1,X));
end
