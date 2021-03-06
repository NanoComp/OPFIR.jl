function Tv(p, N0, N3)
    if ~(N0/N3 >0)
      p.err_tv = true
    #   println(p.alpha_r)
      return 0.
    end
    T = -(p.EG-p.E3) / (p.kB*log(N0/N3)) /8065.73
    return T
end
function fraction_Vg(p, T)
    Q = Qv_0(p, T)
    return 1.0/Q
end

function fraction_V3(p, T)
    Q = Qv_0(p, T)
    return exp(-p.E3/(p.kB*T*8065.73))/Q
end

function Qv_0(p, T)
    if p.script == 1
      data = readdlm("/Users/fanwang/.julia/v0.4/OPFIR/src/E_vib.jl")
    else
      data = readdlm("/Users/fanwang/.julia/v0.4/OPFIR/src/E_vib.jl")
    end
    Q = 1.0
    for i in 1:size(data, 1)
        Q += data[i,2] * exp(-data[i, 1]/(p.kB*T*8065.73))
    end
    return Q
end

function solve_Tv(Q_v, p)
    Qv_diff(T) = Qv_0(p, T) - Q_v
    function Qv(x, fvec)
        fvec[1] = Qv_diff(x[1])
    end
    T_v = nlsolve(Qv, [300.0])
    return T_v.zero
end

function updateTv(p, sol)
    for j in 1:p.num_layers
        N0A = (sol[p.layer_unknown*j-5] + p.ntotal*p.f_G_0/2)
        N3A = (sol[p.layer_unknown*j-4] + p.ntotal*p.f_3_0/2)
        # for vi in 1:p.num_freq
        #   index_b = p.layer_unknown*(j-1) + 1 + (vi-1)*p.n_rot
        #   index_e = index_b + p.n_rot÷2 - 1
        #   N0A += sum(sol[index_b:index_e])
        #   index_b += p.n_rot÷2
        #   index_e += p.n_rot÷2
        #   N3A += sum(sol[index_b:index_e])
        # end
        p.T_vA[j] = Tv(p, N0A, N3A)

        N0E = (sol[p.layer_unknown*j-2] + p.ntotal*p.f_G_0/2)
        N3E = (sol[p.layer_unknown*j-1] + p.ntotal*p.f_3_0/2)
        p.T_vE[j] = Tv(p, N0E, N3E)

        if p.err_tv == true
          println("Tv error is detected! Return current values")
          println("current absorption coefficient: ", p.alpha_r)
          println("problem comes from layer: ", j)
          println("N0A=", N0A, ", N3A=", N3A, ", N0E=", N0E, ", N3E=", N3E)
        #   exit()
          return 0
        end
    end
    # N0A = 0
    # N3A = 0
    # for j in 1:p.num_layers
    #     N0A += (sol[p.layer_unknown*j-5] + p.ntotal*p.f_G_0/2) * p.r_int[j]
    #     N3A += (sol[p.layer_unknown*j-4] + p.ntotal*p.f_3_0/2) * p.r_int[j]
    # end
    # T_v = Tv(p, N0A, N3A)
    # p.T_vA[:] = T_v
    # p.T_vE[:] = T_v
    # println(T_v)
end

function updateks(p)
    for j in 1:p.num_layers
        p.f_GA[j] = fraction_Vg(p, p.T_vA[j])
        p.f_3A[j] = fraction_V3(p, p.T_vA[j])
        p.f_6A[j] = 1 - p.f_GA[j] - p.f_3A[j]
        p.k36A[j] = 1000.0
        p.k63A[j] = p.k36A[j] * p.f_3A[j]/p.f_6A[j]
    end

    for j in 1:p.num_layers
        p.f_GE[j] = fraction_Vg(p, p.T_vE[j])
        p.f_3E[j] = fraction_V3(p, p.T_vE[j])
        p.f_6E[j] = 1 - p.f_GE[j] - p.f_3E[j]
        p.k36E[j] = 1000.0
        p.k63E[j] = p.k36E[j] * p.f_3E[j]/p.f_6E[j]
    end
end
