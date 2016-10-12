using JLD
using ODE

function func(p; sol_start=Array[])
    # initiate some of the parameters from alpha_0
    if p.solstart_flag==0
        sol_0 = zeros(p.num_layers * p.layer_unknown)
    else
        p_sol = load("./p_sol.jld")
        sol_0 = p_sol["sol"]
        tmp_power = p.power
        p = p_sol["p"]
        p.power = tmp_power
        #sol_0 = sol_start
    end

    println("matrix size is ", size(sol_0,1))

    T_err = 1.0
    alpha_err = 1.0
    # sol_0 = zeros(p.num_layers * p.layer_unknown)

    # rel_err = Float64[]
    # sol_0 = andersonaccel(x -> begin
    #         y = fixedpoint(x, p)
    #         push!(rel_err, norm(y - x) / norm(y))
    #         y
    #       end, sol_0, reltol=1e-4)

    while T_err > 5e-2 && alpha_err > 1e-3
        # if p.solstart_flag==0
        #     sol_0 = zeros(p.num_layers * p.layer_unknown)
        # else
        #     sp_sol = load("./p_sol.jld")
        #     sol_0 = p_sol["sol"]
        # end
        alpha_p = deepcopy(p.alpha_r)
        rel_err = Float64[]
        sol_0 = andersonaccel(x -> begin
                y = fixedpoint(x, p)
                push!(rel_err, norm(y - x) / norm(y))
                y
              end, sol_0, reltol=1e-4, m=40)
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
        println("T error: ", T_err, ", alpha_err:", alpha_err)
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
