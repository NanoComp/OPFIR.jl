using Pkg
Pkg.activate("../")
import OPFIR

function paramsweep()
    THzpower = []
    pressures = 10:2:60
    for pressure in pressures
        p = OPFIR.N2O(pressure=pressure,
        power=0.25, num_layers=6, JL=14, n_vib=10,σ_DD=35, σ_GKC=15, σ_SPT=15,qclbroadening=1.e6)
        a, solnew, pnew, taus = OPFIR.outputpower(p, "U", "TE01")
        push!(THzpower, a)
    end
    return pressures, THzpower
end

paramsweep()
