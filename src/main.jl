function func(pressure, power, num_layers)

    p = Params(pressure=pressure, power=power, num_layers=num_layers)
    sol_0 = zeros(p.num_layers * p.layer_unknown)

    println("matrix size is ", size(sol_0,1))
    rel_err = Float64[]
    sol_0 = andersonaccel(x -> begin
            y = fixedpoint(x, p)
            push!(rel_err, norm(y - x) / norm(y))
            y
        end, sol_0, reltol=1e-6)
    return (p, sol_0, rel_err)
end
