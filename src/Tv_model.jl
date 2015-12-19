function Tv(p, N0, N3)
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
    data = readdlm("../src/E_vib.jl")
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
