function Tv(p, N0, N3)
    T = (p.EG-p.E3) / (p.kB*log(N0/N3)) /8065.73
    return T
end

function Qv(kB, T)
    data = readdlm("../src/E_vib.jl")
    Q = 1.0
    for i in 1:size(data, 1)
        Q += data[i,2] * exp(-data[i, 1]/(kB*T*8065.73))
    end
    return Q
end

function fraction_Vg(p, T)
    Q = Qv(p.kB, T)
    return 1.0/Q
end

function fraction_V3(p, T)
    Q = Qv(p.kB, T)
    return exp(-p.E3/(p.kB*T*8065.73))/Q
end
