function func(p; sol_start=Array[])
    # p = Params(pressure=pressure, power=power, num_layers=num_layers, model_flag=model_flag)
    if p.solstart_flag==0
        sol_0 = zeros(p.num_layers * p.layer_unknown)
    else
        sol_0 = sol_start
    end

# update the value of k63 and k36
    # if p.solstart_flag==0
    # else
    #     tot_pump = 0.0
    #     T_vA = Tv(p, sol_0[end-p.n_vib+1]+p.ntotal*p.f_G_0/2,
    #                  sol_0[end-p.n_vib+2]+p.ntotal*p.f_3_0/2)
    #     Q_v = Qv_0(p, T_vA)
    #     println("Tv = ", T_vA)
    #     p.f_G[:] = fraction_Vg(p, T_vA)
    #     p.f_3[:] = fraction_V3(p, T_vA)
    #     p.f_6[:] = 1 - p.f_G[:] - p.f_3[:]
    #     for ri in 1:p.num_layers
    #         # pump_ri = pump_total(p, sol_0, ri)*p.f_6[ri]/(p.f_3[ri]+p.f_6[ri])
    #         # Q_v = p.Q/(1-pump_total(p, sol_0, ri)*Q_v/p.kwall[ri]/p.ntotal)
    #         # if kwall very small, Q_v -> negative, not possible
    #         # T_v = solve_Tv(Q_v, p)
    #         p.netrate_36A[ri] = pump_total(p, sol_0, ri) * p.f_6[ri]*Q_v/(Q_v-1)
    #     end
    # end
    # 
    println("matrix size is ", size(sol_0,1))

    rel_err = Float64[]
    sol_0 = andersonaccel(x -> begin
            y = fixedpoint(x, p)
            push!(rel_err, norm(y - x) / norm(y))
            y
        end, sol_0, reltol=1e-6)
    return (p, sol_0, rel_err)
end
