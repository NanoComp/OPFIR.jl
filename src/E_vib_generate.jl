using DelimitedFiles

function E_vib_generate(name)
    if name == "HCN"
        E1, E2, E3 = 3311.0, 712.0, 2097.0
    elseif name == "N2O"
        E1, E2, E3 = 2224.0, 589.0, 1285.0
    end
    g1, g2, g3 = 1,2,1
    filepath = joinpath(@__DIR__, "viblevels_"*name*".txt")
    if isfile(filepath)
        return
    end
    Es = []
    gs = []
    for i = 0:100
        for j = 0:100
            for k = 0:100
                Eijk = i*E1 + j*E2 + k*E3
                if Eijk < 200*max(E1, E2, E3)
                    push!(Es, Eijk)
                    j>0 ? push!(gs, 2) : push!(gs, 1)
                end
            end
        end
    end
    p = sortperm(Es)
    open(filepath, "w") do io
           writedlm(io, [Es[p][2:end] gs[p][2:end]], ' ')
       end

end
