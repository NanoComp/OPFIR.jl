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
    data = viblevelsN2O()
    Q = 1.0
    for i in 1:size(data, 1)
        Q += data[i,2] * exp(-data[i, 1]/(p.kB*T*8065.73))
    end
    return Q
end

function updateTv(p, sol)
    for j in 1:p.num_layers
        N0 = totN0r(p, sol, j)
        N3 = totN3r(p, sol, j)
        # p.T_vA[j] = Tv(p, N0, N3)
        p.T_vA[j] = p.T
        if p.err_tv == true
          println("Tv error is detected! Return current values")
          println("current absorption coefficient: ", p.alpha_r)
          println("problem comes from layer: ", j)
          println("N0=", N0, ", N3=", N3)
          return 0
        end
    end
end

function updateks(p)
    for j in 1:p.num_layers
        p.f_GA[j] = fraction_Vg(p, p.T_vA[j])
        p.f_3A[j] = fraction_V3(p, p.T_vA[j])
        p.f_6A[j] = 1 - p.f_GA[j] - p.f_3A[j]
        p.k36A[j] = 1000.0
        p.k63A[j] = p.k36A[j] * p.f_3A[j]/p.f_6A[j]
    end

end
