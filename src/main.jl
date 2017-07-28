using JLD
using ODE

function func(p; sol_start=Array[])
    # initiate some of the parameters from alpha_0
    if p.solstart_flag==0
        sol_0 = zeros(p.num_layers * p.layer_unknown + p.n_vib)
        p_0 = p
        matrix_0 = 0
        lu_mat0 = 0
    else
        tmp_power = p.power
        p_sol = load("./p_sol.jld")
        sol_0 = p_sol["sol"]
        sol_in = deepcopy(sol_0)
        p_0 = p_sol["p"]
        p = deepcopy(p_0)
        p.power = tmp_power
        p.solstart_flag = 1
    end

    println("matrix size is ", size(sol_0,1))

    T_err = 1.0
    N_err = 1.0
    while T_err > 2e-2 && N_err > 1e-3
        N_prev = total_NuInv(p, sol_0)
        alpha_p = deepcopy(p.alpha_r)
        rel_err = Float64[]
        sol_in = deepcopy(sol_0)

        if p.solstart_flag == 1
            max_ele = p.num_freq * p.num_layers * (p.n_rot*(p.n_rot+2) + p.n_vib*(p.n_rot+p.n_vib+2))
            rowind_0 = ones(Int64, max_ele)
            colind_0 = ones(Int64, max_ele)
            value_0 = zeros(max_ele)
            compute_row_col_val(rowind_0, colind_0, value_0, p_0, sol_in)
            matrix_0 = sparse(rowind_0, colind_0, value_0)
            mat_modify(matrix_0, p)
            lu_mat0 = lufact(matrix_0)
        end

        sol_0 = andersonaccel(x -> begin
                y = fixedpoint(x, p, matrix_0, lu_mat0)
                push!(rel_err, norm(y - x) / norm(y))
                y
            end, sol_0, reltol=1e-6, m=40)

        if p.model_flag == 1
           break
        end

        alpha_err = norm(alpha_p - p.alpha_r)/norm(p.alpha_r)

        T1 = [deepcopy(p.T_vA); deepcopy(p.T_vE)]
        updateTv(p, sol_0)
        if p.err_tv == true
          return (p, sol_0)
        end
        updateks(p)
        T2 = [deepcopy(p.T_vA); deepcopy(p.T_vE)]
        T_err = norm(T1-T2)/norm(T2)

        N_curr = total_NuInv(p, sol_0)
        N_err = abs(N_curr - N_prev)/abs(N_curr)
        println("T error: ", T_err, ", N_err:", N_err)
        flush(STDOUT)
    end

    return (p, sol_0)
end


function func_tevol(p)
  sol_0 = zeros(p.num_layers * p.layer_unknown)
  sol_f = zeros(p.num_layers * p.layer_unknown, length(p.evol_t)-1)

  dndt(t, sol) = dndt_p(p, t, sol)

  for i in 1:length(p.evol_t)-1
    println("t = ", p.evol_t[i])
    println("L_eff = ", p.L_eff)
    println("alpha = ", sum(p.alpha_r.*p.r_int)/sum(p.r_int))
    tspan = [p.evol_t[i]; p.evol_t[i+1]]
    (t, sol_total) = ode45(dndt, sol_0, tspan; abstol=1e-3, reltol=1e-2)
    sol_0 = sol_total[end]
    sol_f[:, i] = sol_total[end]

    if p.model_flag == 2
      updateTv(p, sol_0)
      if p.err_tv == true
        break
      end
      updateks(p)
    end
  end

  return (p, sol_f)
end


function dndt_p(p, t, sol)
  max_ele = p.num_freq * p.num_layers * (p.n_rot*(p.n_rot+2) + p.n_vib*(p.n_rot+p.n_vib+2))

  rowind = ones(Int64, max_ele)
  colind = ones(Int64, max_ele)
  value = zeros(max_ele)
  rhs = zeros(p.num_layers*p.layer_unknown)

  update_alpha_from_N!(p, sol)
  # println(p.alpha_r)
  update_Param_from_alpha!(p, sol)

  compute_rhs(rhs, p, sol)
  compute_row_col_val(rowind, colind, value, p, sol)

  matrix = sparse(rowind, colind, value)

  return (matrix*sol + rhs)

end

function total_NuInv(p, sol)
    invU = zeros(p.num_layers)
    for k = 1:p.num_layers
        invU[k] = sum(inv_U_dist_layer(p, sol, k))
        # invU[k] = Nu_NT_dist_layer(p, sol, k)[1] - p.g_U/p.g_L*Nu_1_NT_dist_layer(p, sol, k)[1]
    end
    return sum(p.r_int.*invU)/sum(p.r_int)
end

function totinv(p, sol, llevel)
    popinv = zeros(p.num_layers)
    if llevel == 'U'
        for k = 1:p.num_layers
            popinv[k] = sum(inv_U_dist_layer(p, sol, k))
        end
    elseif llevel == 'L'
        for k = 1:p.num_layers
            popinv[k] = sum(inv_L_dist_layer(p, sol, k))
        end
    else
        throw(ArgumentError("Lasing level error!"))
    end
    return sum(p.r_int.*popinv)/sum(p.r_int)
end


function zerobessel(m)
    # mode_num: 1: TE01 / 2: TE12 / 3: TE02 / 4: TE22 / 5: TE11 / 6: TE21 / 7: TM01 / 8: TM11
    #zeros of Bessel functions:
    # p_library = [3.83, 5.33, 7.02, 6.71, 1.84, 3.05, 2.4, 3.83]
    return m == "TE01" ? 3.83 :
           m == "TE12" ? 5.33 :
           m == "TE02" ? 7.02 :
           m == "TE22" ? 6.71 :
           m == "TE11" ? 1.84 :
           m == "TE21" ? 3.05 :
           m == "TM01" ? 2.40 :
           m == "TM11" ? 3.83 :
           throw(ArgumentError("Cavity mode entered is not supported yet!"))
end

function totinvmode(p, sol, llevel, cavitymode)
    m = parse(Int, cavitymode[3])
    radius_m = p.radius/100
    kρ = zerobessel(cavitymode)/radius_m
    ϕ_list = linspace(0, 2*pi, 150)
    denom = 0.
    numerator = 0.
    pop_inv = zeros(p.num_layers)

    for ri = 1:p.num_layers
        pop_inv[ri] = llevel=='U' ? sum(inv_U_dist_layer(p, sol, ri)) :
                      llevel=='L' ? sum(inv_L_dist_layer(p, sol, ri)) :
                      throw(ArgumentError("Lasing level error!"))
    end
    E_sq = 0.
    for ri = 1:p.num_layers
        r = p.r_int[ri]
        for ϕ in ϕ_list
            if contains(cavitymode, "TM") ## TM modes
                Ez = besselj(m, kρ*r) * cos(m*ϕ)
                Er = derv_bessel(m, kρ*r) * cos(m*ϕ)
                Eϕ = m/(kρ*r) * besselj(m, kρ*r) * (-sin(m*ϕ))
            elseif contains(cavitymode, "TE") ## TE modes
                Ez = 0
                Er = m/(kρ*r) * besselj(m, kρ*r) * (-sin(m*ϕ))
                Eϕ = derv_bessel(m, kρ*r) * cos(m*ϕ)
            else
                throw(ArgumentError("cavity mode can only be TE or TM"))
            end
            E_sq = Ez*conj(Ez) + Er*conj(Er) + Eϕ*conj(Eϕ)

            # E_sq = 1.
            denom += E_sq * r # cylindrical coordinate
            numerator += E_sq * r * pop_inv[ri]
        end
        # println(E_sq)
    end

    return numerator/denom
end

function total_NlInv(p, sol)
    invL = zeros(p.num_layers)
    for k = 1:p.num_layers
        invL[k] = sum(inv_L_dist_layer(p, sol, k))
    end
    return sum(p.r_int.*invL)/sum(p.r_int)
end

function nonthpopinv(p, sol)
    invU = zeros(p.num_layers)
    invL = zeros(p.num_layers)
    for kl = 1:p.num_layers
        invU[kl] = sum(Nu_NT_dist_layer(p,sol,kl) - p.g_U/p.g_L*Nu_1_NT_dist_layer(p,sol,kl))
        invL[kl] = sum(Nl_1_NT_dist_layer(p,sol,kl) - p.g_U/p.g_L*Nl_NT_dist_layer(p,sol,kl))
    end
    total_invU = sum(p.r_int.*invU)/sum(p.r_int)
    total_invL = sum(p.r_int.*invL)/sum(p.r_int)
    return (total_invU, total_invL)
end

function popinvth(p)
    ν0 = (p.f_dir_lasing+p.f_ref_lasing)/2
    lambda = p.c/ν0
    alpha = cavityloss(p) # unit m^-1
    Nt = 8*pi^2 / lambda^2 * p.t_spont * alpha * sum(p.Δ_f_NTF.*p.r_int)/sum(p.r_int)
    return Nt
end

function ohmicloss(p, llevel, cavitymode)
# p_library = [2.40 3.83 1.84 3.05 3.83]
# % 1 -> TM01 cutoff
# % 2 -> TM11 cutoff
# % 3 -> TE11 cutoff
# % 4 -> TE21 cutoff
# % 5 -> TE01 cutoff

# mode_num: 1: TE01 / 2: TE12 / 3: TE02 / 4: TE22 / 5: TE11 / 6: TE21 / 7: TM01 / 8: TM11
#zeros of Bessel functions:
# p_library = [3.83, 5.33, 7.02, 6.71, 1.84, 3.05, 2.4, 3.83]
radius_m = p.radius/100
# x0 = (llevel=='U' ? p.p_library[3] : p.p_library[1])# TE01 mode
x0 = zerobessel(cavitymode)
m = parse(Int, cavitymode[3])

resitivityCu = 1/4e7 # copper resistivity at room temperature;
conductivityCu = 1/resitivityCu # copper conductivity at room temperature;
mu = 4*pi*1e-7 # magnetic permeability in copper

eta = 377 # impedance of air in ohms
f0 = (llevel=='U' ? p.f_dir_lasing : p.f_ref_lasing)
lambda = p.c/f0 # wavelength in m
k0 = 2*pi/lambda # wavevector
Rs = sqrt(pi*f0*mu/conductivityCu) # in ohms

if contains(cavitymode, "TE")
    lossval = Rs/(radius_m*eta*sqrt(1-(x0/k0/radius_m)^2))*
            ((x0/k0/radius_m)^2+m^2/(x0^2-m^2)) #; % in 1/m
elseif contains(cavitymode, "TM")
    lossval = Rs/(radius_m*eta*sqrt(1-(x0/k0/radius_m)^2))
end
return lossval
end


function cavityloss(p, llevel, cavitymode)
    Rback = 1.
    Rfront = 0.96 * Rback
    alpha = 2 * ohmicloss(p, llevel, cavitymode) - log(Rfront*Rback)/(2p.L/100)
    return alpha
end


function gaincoeffcient(f, Φ, p, sol, taus, level)
    γ = 0.0
    for vi in 1:p.num_freq
        if level == 'L'
            f0 = p.f_dist_ref_lasing[vi]
        elseif level == 'U'
            f0 = p.f_dist_dir_lasing[vi]
        else
            throw(ArgumentError("level can only be L or U!"))
        end
            # popinv: spatial averged pop inversion for U
        popinv = 0.
        if level == 'L'
            for i = 1:p.num_layers
                popinv += (OPFIR.inv_L_dist_layer(p, sol, i)[vi]*p.r_int[i])/sum(p.r_int)
            end
        else
            for i = 1:p.num_layers
                popinv += (OPFIR.inv_U_dist_layer(p, sol, i)[vi]*p.r_int[i])/sum(p.r_int)
            end
        end
        Δnu = p.Δ_fP
        λ = p.c/f
        sigma_f = λ^2/8/π/p.t_spont * 1/pi * Δnu/((f-f0)^2 + Δnu^2)
        if Φ == 0
            γ += popinv * sigma_f
        else
            Φs = 1/taus/sigma_f
            γ += popinv * sigma_f / (1+Φ/Φs)
        end
    end
    return γ
end

function comptaus(nonth_popinv, wi_list)
    taus = similar(nonth_popinv[:, :, 1])
    for pressure in 1:size(nonth_popinv, 1)
        for power in 1:size(nonth_popinv, 2)
            tmp = vec(nonth_popinv[pressure, power, :]) / nonth_popinv[pressure, power, 1]
            a, b = linreg(wi_list*1.0, 1./tmp - 1)
            taus[pressure, power] = b
        end
    end
    return taus
end


function compfraction(p, sol)
    N0A = N3A = NΣA = 0.
    for j in 1:p.num_layers
        N0A += ((sol[p.layer_unknown*j-5] + p.ntotal*p.f_G_0/2)) * p.r_int[j]
        N3A += ((sol[p.layer_unknown*j-4] + p.ntotal*p.f_3_0/2)) * p.r_int[j]
        NΣA += ((sol[p.layer_unknown*j-3] + p.ntotal*p.f_6_0/2)) * p.r_int[j]
    end
    f0 = N0A/sum(N0A+N3A+NΣA)
    f3 = N3A/sum(N0A+N3A+NΣA)
    fΣ = NΣA/sum(N0A+N3A+NΣA)
    return (f0, f3, fΣ)
end
