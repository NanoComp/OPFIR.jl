function func(p; sol_initial=Array[], initial_flag=0)
    # p = Params(pressure=pressure, power=power, num_layers=num_layers, model_flag=model_flag)
    if initial_flag==0
        sol_0 = zeros(p.num_layers * p.layer_unknown)
    else
        sol_0 = sol_initial
    end

    if initial_flag==0
    else
        tot_pump = 0.0
        temp_36 = 0.0
        for ri in 1:p.num_layers
            tot_pump += pump_total(p, sol_0, ri)*p.f_6_0/(p.f_3_0+p.f_6_0)*ri
            ind3 = p.layer_unknown*ri - p.n_vib + 2
            temp_36 += (sol_0[ind3] - p.f_3_0/p.f_6_0 * sol_0[ind3+1])*ri
        end
        p.k63 = tot_pump/temp_36
        p.k36 = p.f_6_0/p.f_3_0 * p.k63
    end
    println("k36 = ", p.k36, " k63 = ", p.k63)
    println("matrix size is ", size(sol_0,1))

    rel_err = Float64[]
    sol_0 = andersonaccel(x -> begin
            y = fixedpoint(x, p)
            push!(rel_err, norm(y - x) / norm(y))
            y
        end, sol_0, reltol=1e-6)
    return (p, sol_0, rel_err)
end
