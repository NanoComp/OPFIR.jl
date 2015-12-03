function func(pressure, power, num_layers, niter::Int64)
    ######################## calculate parameters ##########################
    para = OPFIR_compute_parameters(pressure, power, num_layers);
    layer_unknown = n_rot*para.num_freq + n_vib
    para.niter = niter;
#     para.D = 0;
    ######################### create matrix and rhs ###################

    sol_0 = zeros(para.num_layers * layer_unknown);
    println("matrix size is ", size(sol_0,1))
    rel_err = Float64[]
    sol_0 = andersonaccel(x -> begin
            y = OPFIR_fixedpoint(x, para)
            push!(rel_err, norm(y - x) / norm(y))
            y
        end, sol_0, reltol=1e-6)
    return (para, sol_0, rel_err);
end
