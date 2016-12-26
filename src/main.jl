using JLD
using ODE

function func(p; sol_start=Array[])
    # initiate some of the parameters from alpha_0
    if p.solstart_flag==0
        sol_0 = zeros(p.num_layers * p.layer_unknown)
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
    end
    return sum(p.r_int.*invU)/sum(p.r_int)
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

function cavityloss(p)
    Q = 10000
    ν0 = (p.f_dir_lasing+p.f_ref_lasing)/2
    lambda = p.c/ν0
    L = p.L/100 # meter
    t_rt = 2*L/p.c # in sec
    alpha_perc = 2*pi*ν0 * t_rt/Q

    alpha = -log(1-alpha_perc)/2/L # unit m^-1
    return alpha
end
