using DataFrames, Polynomials

function getoutpower(dirlabel;
    plist = vcat(5:5:200),
        powerlist = vcat(0.27),
        Jlist = vcat(16),
        DDlist = vcat(15.:10:85),
        GKClist = vcat(15:10.:85),
        SPTlist = vcat(15:10:85.),
        nlayers = 6,
        lossfactor = 1.0)
    """
        returns a dataframe:
        pressure, pump power, THz power, DD, GKC, SPT, J, spontaneous emission time, Δν_P, gain coefficient
    """
    cd(homedir()*"/.julia/v0.6/OPFIR")
    cd(string(dirlabel))
    outpower_DD_GKC = DataFrame(pressure=Float64[],
        pumppower=Float64[],
        THzpower=Float64[],
        DD=Float64[],
        GK=Float64[],
        SPT=Float64[],
        J=Float64[],
        spont = Float64[],
        Δν_P = Float64[],
        gain = Float64[],
        lossfactor = Float64[])

    for iDD in 1:length(DDlist)
        DD = DDlist[iDD]
        for iGKC in 1:length(GKClist)
            GKC = GKClist[iGKC]
            for iSPT in 1:length(SPTlist)
                SPT = SPTlist[iSPT]
                THzpower = Array{Float64}(0)
                pressures = similar(THzpower)
                absorption = similar(THzpower)
                cavitylen = similar(THzpower)
                p = 0
                for J in Jlist
                    for power in powerlist
                        for pressure in plist
                            filepath = string("p_",pressure,"_power_",power,"_J_",J,"_nlayers_",nlayers,
                                "_DD_",DD,"_GKC_",GKC, "_SPT_", SPT)
                            if !isdir(filepath)
                               continue
                            end
                            if !isfile(filepath*"/power_results.txt")
                               continue
                            end
                            taus = 0.0
                            x = 0.0
                            open(filepath*"/power_results.txt") do f
                                readline(f)
                                x = float(readline(f))
                                push!(THzpower, x)
                                readline(f)
                                taus = float(readline(f)[6:end])
                            end
                            p_sol = load(filepath*"/p_sol.jld")
                            p = p_sol["p"]
                            sol = p_sol["sol"]
                            # Acoeff = 64*pi^4/3/6.63e-27/(3e10)^3 * p.f_dir_lasing^3 * 0.17^2 * p.JU/(2p.JU+1) * 1e-36
                            # p.t_spont = 1/Acoeff
                            if lossfactor != 1.
                                x = OPFIR.outpowermode(p, sol, "U", "TE01", taus, lossfactor=lossfactor)[1]
                            end
                            gain = OPFIR.gaincoefmode(0, p.f_dir_lasing, p, sol, "U", "TE01", taus)
                            push!(outpower_DD_GKC, (pressure, power, x, DD, GKC, SPT, J, p.t_spont, p.Δ_fP, gain, lossfactor))
                        end
                    end
                end
            end
        end
    end
    return outpower_DD_GKC
end


function getthreshold(irpower, THzpower)
    """ return pump threshold and slope of THz power vs IR power """
    x = Array{Float64}(0)
    y = Array{Float64}(0)
    for i in 1:length(THzpower)
        if THzpower[i] > 0
            push!(x, irpower[i])
            push!(y, THzpower[i])
        end
    end
    if length(y) < 2
        return -1, -1
    end
    if length(y) == 2
        return (y[2]-y[1])/(x[2]-x[1]), -y[1]/((y[2]-y[1])/(x[2]-x[1]))+x[1]
    end
    (y[2]-y[1])/(x[2]-x[1]), Polynomials.roots(polyfit(x[1:3], y[1:3], 2))[2]
end
