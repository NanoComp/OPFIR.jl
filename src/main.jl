function func(p; sol_start=Array[])
    # p = Params(pressure=pressure, power=power, num_layers=num_layers, model_flag=model_flag)
    if p.solstart_flag==0
        sol_0 = zeros(p.num_layers * p.layer_unknown)
    else
        sol_0 = sol_start
    end

    println("matrix size is ", size(sol_0,1))


    T_err = 1.0
    while T_err > 1e-3
        rel_err = Float64[]
        sol_0 = andersonaccel(x -> begin
                y = fixedpoint(x, p)
                push!(rel_err, norm(y - x) / norm(y))
                y
                end, sol_0, reltol=1e-6)

        T1 = [deepcopy(p.T_vA); deepcopy(p.T_vE)]
        updateTv(p, sol_0)
        updateks(p)
        T2 = [deepcopy(p.T_vA); deepcopy(p.T_vE)]
        T_err = norm(T1-T2)/norm(T2)
        println(T_err)
        
        sol_0 = zeros(p.num_layers * p.layer_unknown)
    end

    return (p, sol_0)
end
