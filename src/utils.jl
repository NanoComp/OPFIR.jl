using DataFrames, Polynomials

function readpsol(dirlabel, DD::Float64, GKC::Float64, SPT::Float64, J::Int64, pressure::Int64, power::Float64, nlayer::Int64, llevel)
    """readpsol(dirlabel, DD::Float64, GKC::Float64, SPT::Float64, J::Int64, pressure::Int64, power::Float64, nlayer::Int64, llevel)
    """
    cd(homedir()*"/.julia/v0.6/OPFIR")
    cd(string(dirlabel))
    filepath = string("p_",pressure,"_power_",power,"_J_",J,"_nlayers_",nlayer, "_DD_",DD,"_GKC_",GKC, "_SPT_", SPT)
    if llevel in ['L', "L"]
        filepath = string(filepath, "_L")
    end
    if !isdir(filepath)
        return 0, 0, 0, 0
        throw(ArgumentError("check the path!"))
        return
    end
    if !isfile(filepath*"/power_results.txt")
       throw(ArgumentError("no file exits!"))
    end
    p_sol = load(filepath*"/p_sol.jld")
    p = p_sol["p"]
    sol = p_sol["sol"]

    taus = 0.0
    THzpower = 0.
    open(filepath*"/power_results.txt") do f
        readline(f)
        x = float(readline(f))
        THzpower = x
        readline(f)
        taus = float(readline(f)[6:end])
    end

    return p, sol, THzpower, taus
end


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

    J0 = []
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
                            if !(J in J0)
                                println("J: ", J)
                                println("ohmic loss: " * string(ohmicloss(p, "U", "TE01")) * " m-1")
                                println("transmission/reflection loss: " *
                                string(cavityloss(p, "U", "TE01")-ohmicloss(p, "U", "TE01")) * " m-1")
                                push!(J0, J)
                            end
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

function contourplot(df, J; plotth=false)
    k = 1
    uniquedf = unique(df, [:DD, :SPT, :GK, :J]) # a DataFrame
    nplots = size(uniquedf, 1)
    figure(figsize=(12, 4nplots))
    for udf in eachrow(uniquedf)
        dfi = df[(df[:DD].==udf[:DD]) .&
                 (df[:SPT].==udf[:SPT]) .&
                 (df[:GK].==udf[:GK]) .&
                 (df[:J].==udf[:J]), :]
        pressures = unique(dfi[:pressure])
        pumppowers = unique(dfi[:pumppower])
        THzpower = zeros(length(pressures), length(pumppowers))
        X = similar(THzpower)
        Y = similar(THzpower)
        th = zeros(size(pressures))

        for i in 1:length(pressures)
            for j in 1:length(pumppowers)
                THzpower[i,j] = dfi[(dfi[:pressure].==pressures[i]) .&
                                     (dfi[:pumppower].==pumppowers[j]), :THzpower][1]
                X[i,j] = pressures[i]
                Y[i,j] = pumppowers[j]
                THzpower[i,j] = THzpower[i,j] < 0 ? 0 : THzpower[i,j]
            end
            if plotth
                th[i] = getthreshold(dfi[dfi[:pressure].==pressures[i], :pumppower],
                                     dfi[dfi[:pressure].==pressures[i], :THzpower])[2]
            end
        end
#
#     DDlist = unique(df[:DD])
#     for DD in DDlist
#         pressurelist = vcat(20:10:80.)
#         powerlist = vcat(0.01:0.03:0.25)
#         THzpower = zeros(length(pressurelist), length(powerlist))
#         X = zeros(length(pressurelist), length(powerlist))
#         Y = zeros(length(pressurelist), length(powerlist))
#
#         dfDD = df[(df[:DD].==DD) .& (df[:J].==J), :]
#         for i in 1:length(pressurelist)
#             for j in 1:length(powerlist)
#                 if i > 2
#                     THzpower[i,j] = dfDD[(dfDD[:pressure].==pressurelist[i]) .&
#                                          (dfDD[:pumppower].==powerlist[j]), :THzpower][1]
#                 end
#                 X[i,j] = pressurelist[i]
#                 Y[i,j] = powerlist[j]
#                 THzpower[i,j] = THzpower[i,j] < 0 ? 0 : THzpower[i,j]
#             end
#         end
# #         THzpower[2, :] = -THzpower[4, :] + 2THzpower[3,:]
# #         THzpower[1, :] = -THzpower[3, :] + 2THzpower[2,:]
#         for i in 1:length(pressurelist)
#             nonzeroTHzpower = Array{Float64}(0)
#             nonzeropumppower = similar(nonzeroTHzpower)
#             for j in 1:length(powerlist)
#                 if THzpower[i,j] > 0
#                     push!(nonzeroTHzpower, THzpower[i,j])
#                     push!(nonzeropumppower, powerlist[j])
#                 end
#                 THzpower[i,j] = THzpower[i,j] < 0 ? 0 : THzpower[i,j]
#             end
#             if length(nonzeroTHzpower) > 2
#                 spl = Spline1D(nonzeropumppower, nonzeroTHzpower, k=2)
#                 THzpower[i,:] = spl(powerlist)
#             end
#         end
        subplot(nplots÷2+1,2,k)
        contourf(X, Y, 1000THzpower, 100)
        colorbar()
        plotth && plot(pressures, th, "ro", markersize=6)
        plotth && plot(vcat(20:5:85,88), vcat(45,51,59,66,77,92,108,120,141,154,177,202,223,251,260)/300*0.26, "ko")
        xlabel("pressure")
        ylabel("pump power")
        title("DD: "*string(dfi[:DD][1])*" "*" A2, SPT: " * string(dfi[:SPT][1]) *"A2,  loss factor: " *
        string(dfi[:lossfactor][1]))
        xlim(pressures[1], pressures[end])
        ylim(0.01,0.24)
        k += 1
    end
end
