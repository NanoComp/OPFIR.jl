function LinearMolecule(DefaultT=Float64;
    radius = 0.25, # in cm
    pump_radius = radius, # in cm
    L = 15, # in cm
    L_eff = L,
    T = 300,

    cavitywall = "Cu",
    cavitymode = "TE01",
    frontmirrorT_THz = 0.016,
    frontmirrorT_IR = 0.04,
    backmirrorR_IR = 0.95,
    lossfactor = 1.0,
    lasinglevel = "U",

    ##
    name = "N20",
    σ_GKC = -1, # if given, use the given value, otherwise, compute it
    σ_DD = -1,
    σ_SPT = -1,

    pressure = 20.0,
    ntotal = 9.66e24 * pressure * 1e-3 / T,
    pbroadenwidth = -1,
    pbfactor = 1.0,
    ΔfQCL = 1.e6,

    beta13 = 0.0,

    D_factor = 1.0,

    B = -1,
    DJ = -1,

    JL = 11,
    pumpbranch = "R",
    K0 = 0,

    f_offset = 0,

    dipolemoment = -1,
    dipolederivative = -1,

    nrot_aboveJL = 8,
    n_vib = 10,

    n0 = 1.0,

    power = 0.25,

    num_layers = 6,
    solstart_flag = 0,
    mumps_solver = 0,
    WiU = 0,
    WiL = 0,

    taus = 0.,
    THzpower = 0.,
    gaincoeff = 0.
    )

    ###################################################
    #### physical constants
    ###################################################
    ##
    h = 6.626068e-34
    c = 2.99792458e8
    ev = 1.60217646e-19
    kB = 8.617342e-5 # in ev/K
    norm_time = 1e6 # for unit
    NA = 6.0221413e23
    mu0 = 4e-7*pi
    eps0 = 8.85e-12
    kBT = kB*T*8065.73 # in cm^-1

    ###################################################
    #### molecule
    ###################################################
    if name == "N2O"
        M = 44
        σ_GKC = σ_GKC==-1 ? 15.0 : σ_GKC
        σ_DD = σ_DD==-1 ? 35.0 : σ_DD
        σ_SPT = σ_SPT==-1 ? 15.0 : σ_SPT
        pbroadenwidth = pbroadenwidth==-1 ? 4.0e6 : pbroadenwidth
        B = B==-1 ? 0.419*c/1e7 : B # in GHz
        DJ = DJ==-1 ? 17.6e-8*c/1e7 : DJ
        EG, E3 = 0, 2224.0
        g3 = 1
        dipolemoment = 0.17
        dipolederivative = 0.2457
    elseif name == "HCN"
        M = 27
        pbroadenwidth = pbroadenwidth==-1 ? 18.4e6 : pbroadenwidth
        σ_GKC = σ_GKC==-1 ? 25.0 : σ_GKC
        σ_DD = σ_DD==-1 ? 522.0 : σ_DD
        σ_SPT = σ_SPT==-1 ? 15.0 : σ_SPT
        B = B==-1 ? 44.3159757 : B
        DJ = DJ==-1 ? 87.24e-9 : DJ
        EG, E3 = 0, 1424.0
        g3 = 2
        dipolemoment = 2.98
        dipolederivative = 0.0493
    elseif name == "CO"
        M = 28
        pbroadenwidth = pbroadenwidth==-1 ? 3.2e6 : pbroadenwidth
        σ_GKC = σ_GKC==-1 ? 44.0 : σ_GKC
        σ_DD = σ_DD==-1 ? sqrt(T*M)/3.15*(pbroadenwidth/1e6) : σ_DD
        σ_SPT = σ_SPT==-1 ? 15.0 : σ_SPT
        B = B==-1 ? 57.635968 : B
        DJ = DJ==-1 ? 0.1835058E-6 : DJ
        EG, E3 = 0, 2143.0
        g3 = 1
        dipolemoment = 0.12
        dipolederivative = 0.107555
    end

    Δ_fP0 = pbroadenwidth*(pressure/1000.0) * pbfactor
    Δ_fP = Δ_fP0 + ΔfQCL

    beta13 = 1.20 * sqrt(power)/radius
    Δ_f_RabiF = zeros(num_layers)
    Δ_f_RabiB = zeros(num_layers)
    Δ_f_NTF = ones(num_layers) * Δ_fP
    Δ_f_NTB = ones(num_layers) * Δ_fP

    v_avg = 205*sqrt(T/M) # in m/sec
    MFP = 0.732*T/pressure/σ_GKC # in cm
    # diffusion coefficient in m^2/microsec. Einstein-Smoluchowski equation;
    vel = v_avg/sqrt(2)/norm_time # avg absolute vel in m/microsec
    D = 1/3 * vel * MFP * 1e-2 * D_factor

    if pumpbranch == "R"
        JU = JL+1
    elseif pumpbranch == "P"
        JU = JL-1
    else
        throw(ArgumentError("pump branch can only be P or R!"))
    end
    g_L = 2JL + 1
    g_U = 2JU + 1
    CL = compC(B, DJ, M, T, JL)
    CL1 = compC(B, DJ, M, T, JL+1)
    CU = compC(B, DJ, M, T, JU)
    CU1 = compC(B, DJ, M, T, JU-1)

    EU = E3*c/1e7 + B*JU*(JU+1) - DJ*JU^2*(JU+1)^2 # in GHz
    EL = B*JL*(JL+1) - DJ*JL^2*(JL+1)^2
    f₀ = (EU - EL) * 1e9 # in Hz
    f_pump = f₀ + f_offset
    Δ_f₀D = 3.58e-7*f₀*sqrt(T/M)

    Q = Qv(T, name)
    f_G_0 = exp(-EG/kBT)/Q
    f_3_0 = g3*exp(-E3/kBT)/Q

    # alpha_0 in m^-1
    totalNL0 = CL * ntotal * f_G_0
    totalNU0 = CU * ntotal * f_3_0
    alpha_0 = exp(-log(2)*(f_offset/Δ_f₀D)^2) * sqrt(log(2)/pi)/Δ_f₀D *
              8*pi^3/3/h/c * (totalNL0 - totalNU0) * 1e-36 *
              (dipolederivative^2 * max(JL, JU)/(2JL+1)) * f_pump * 1e-13
    alpha_r = zeros(num_layers)

    f_dir_lasing = (2B*JU - 4DJ*JU^3) * 1e9 # U to U-1
    f_ref_lasing = (2B*(JL+1) - 4DJ*(JL+1)^3) * 1e9 # L+1 to L

    EinsteinA = 64*pi^4/3/(h*1e7)/(c*100)^3 * f_dir_lasing^3 *
                dipolemoment^2 * max(JL, JU)/(2JL+1) * 1e-36
    t_spont = 1/EinsteinA

    Js, n_rot = Jlevels(JL, JU, K0)

    kDD = 19.8 * pressure * σ_DD/ sqrt(T*M)/1000.0 # 1/microsec
    ka = 19.8 * pressure * σ_SPT/ sqrt(T*M)/1000.0
    # in 1/microsec. rate_ij = rate_DD * prob. of ij collision, also kDDmat[i,j]:
    kDDmat = zeros(n_rot, n_rot)
    for i in 2:length(Js)
        kDDmat[i, i-1] = kDDmat[i+n_rot÷2, i+n_rot÷2-1] =
        kDD*Q_selectn_hl(Js[i], K0)
    end
    for i in 1:length(Js)-1
        dE = (2B*(Js[i]+1) - 4DJ*(Js[i]+1)^3) * 1e9 # in Hz
        kDDmat[i, i+1] = kDDmat[i+n_rot÷2, i+n_rot÷2+1] =
        kDD*Q_selectn_lh(Js[i], K0) * exp(-dE * h/1.38e-23/T)
    end

    # cj: vector of fraction of Js in its thermal pool
    ck = zeros(n_rot÷2)
    for k in 1:n_rot÷2
        ck[k] = compC(B, DJ, M, T, Js[k])
    end
    ck = ck/sum(ck)
    cj = vcat(ck, ck)

    vibdata = viblevels(name)
    f_vib = zeros(n_vib)
    f_vib[1] = 1.0/Q
    f_vib[2:end] = vibdata[1:n_vib-1,2] .* exp.(-vibdata[1:n_vib-1, 1]/kBT)/Q
    frel_vib = f_vib/sum(f_vib) # sum to be 1.0; used for BC redistribution
    ## VV collision:
    σ_VV = σ_GKC
    kVV = 19.8 * pressure * σ_VV/ sqrt(T*M) / 1000 # 1/microsec
    kVVmat = zeros(n_vib, n_vib)
    energyV = vcat(0, vibdata[1:n_vib-1, 1])
    gV = vcat(1, vibdata[1:n_vib-1, 2])
    # kVVmat[i,j] denotes transition from i to j, first only consider neighboring v levels
    for i in 2:n_vib
        kVVmat[i, i-1] = kVV * exp(- abs(energyV[i]-energyV[i-1])/kBT)
    end
    for i in 1:n_vib-1
        kVVmat[i, i+1] = kVVmat[i+1, i] * gV[i+1]/gV[i] * exp(- (energyV[i+1]-energyV[i])/kBT)
    end

    ######################################################################
    ### pump power related
    ######################################################################
    powerF = zeros(num_layers)
    powerB = zeros(num_layers)
    averagePF = power * ones(num_layers)
    averagePB = power * ones(num_layers)
    # in unit m^-3 microsec^-1
    pump_0 = 9.4e13 * power/(radius^2)/Δ_f₀D *
             (dipolederivative^2*max(JL, JU)/(2JL+1)) *
             exp(-log(2)*((f_pump-f₀)/Δ_f₀D)^2)/norm_time

    nfP = max(120, 45*150/pressure) # for 45mTorr, 150Δ_fP is enough
    nfP = min(700, nfP) # if pressure is too small, fix as 350Δ_fP
    if pressure <= 30
        nfP = nfP ÷ 2
    end
    f_range =  nfP * Δ_fP
    f_range = min(f_range, 2Δ_f₀D) # if f_range obtained from Δ_fP is too large, use 2Δ_f₀D as half width

    num_freq = round(Int64,max(50,2f_range/(Δ_fP/4))) # resolution is Δ_fP/4
    num_freq = 100

    pump_IR = zeros(num_freq, num_layers)
    layer_unknown = n_rot*num_freq + n_vib

    df = 2.0 * f_range / num_freq
    f_dist_end = range(-f_range, stop=f_range, length=num_freq + 1) .+ f₀
    f_dist_ctr = f_dist_end[1:end-1] .+ df/2

    velocity = (f_dist_ctr .- f₀)/f₀ # in unit c
    f_dist_dir_lasing = velocity * f_dir_lasing .+ f_dir_lasing
    f_dist_ref_lasing = velocity * f_ref_lasing .+ f_ref_lasing
    f_dist_ctrB = f₀ .- f₀ * velocity

    norm_dist = Normal(f₀, Δ_f₀D / sqrt(2*log(2)))
    pdf1 = pdf.(norm_dist, f_dist_ctr)
    gauss_dist = pdf1 * df

    SHBF = zeros(num_freq, num_layers)
    SHBB = zeros(num_freq, num_layers)

    Δr = radius/100 / num_layers # in m
    r_ext = range(0,stop=radius/100, length=num_layers+1)
    r_int = 0.5*(r_ext[1:end-1] + r_ext[2:end]) # in m

    cavityloss = 0

    return ParamsLinearMolecule{DefaultT}(
    radius, pump_radius, L, L_eff, T, cavitywall, cavitymode, frontmirrorT_THz,
    frontmirrorT_IR, backmirrorR_IR, cavityloss, lossfactor, lasinglevel,
    h, c, ev, kB, norm_time, NA, mu0, eps0, kBT,
    name, M, σ_GKC, σ_DD, σ_SPT,
    pressure, ntotal,
    kDD, kDDmat, ka,
    pbroadenwidth, pbfactor, Δ_fP0, ΔfQCL, Δ_fP,
    beta13, Δ_f_RabiF, Δ_f_RabiB, Δ_f_NTF, Δ_f_NTB,
    v_avg, MFP, D, D_factor,
    B, DJ, JL, JU, K0, pumpbranch, CL, CL1, CU, CU1, g_L, g_U,
    f₀, Δ_f₀D, f_offset, f_pump,
    dipolemoment, dipolederivative,
    alpha_0, alpha_r,
    EG, E3, Q, f_G_0, f_3_0,
    f_dir_lasing, f_ref_lasing,
    t_spont,
    Js, n_rot, n_vib, cj, frel_vib, kVVmat, n0,
    power, powerF, powerB, averagePF, averagePB, pump_0, pump_IR,
    f_range, velocity, f_dist_end, f_dist_ctr, f_dist_ctrB, df, num_freq, layer_unknown,
    f_dist_dir_lasing, f_dist_ref_lasing, gauss_dist, SHBF, SHBB,
    Δr, r_int, num_layers,
    solstart_flag, mumps_solver,
    WiU, WiL,
    taus, THzpower, gaincoeff
    )
end

function compC(B, DJ, M, T, JL; nJ=1000)
    # compute fraction of JL in that vibrational level
    Js = 0:1:nJ
    Q = ql = 0.0
    # B, DJ = [0.419, 17.6e-8]*2.99792458e8/1e7 # in GHz
    for J in Js
        Ej = B*J*(J+1) - DJ*J^2*(J+1)^2
        qj = (2J+1) * exp(-Ej*6.626068e-2/(1.38064852*T))
        Q += qj
        if J == JL
            ql = qj
        end
    end
    return ql/Q
end

function Qv(T, name)
    data = viblevels(name)
    Q = 1.0
    for i in 1:size(data, 1)
        Q += data[i,2] * exp(-data[i, 1]/(8.617342e-5*T*8065.73))
    end
    return Q
end

function Jlevels(JL, JU, K0)
    Jmin = 0 #max(K0, JL-N)
    Jmax = JL + 8
    return Jmin:Jmax, (Jmax-Jmin+1)*2
end

function Q_selectn_hl(J, K)
    return (J^2-K^2)/(J*(2*J+1))
end

function Q_selectn_lh(J, K)
    return ((J+1)^2-K^2)/((J+1)*(2*J+1))
end

function f_NT_normalized(ν, Δ_f_NT, f_pump, df)
    SHB = 1/π * df * Δ_f_NT ./ ((ν .- f_pump).^2 .+ Δ_f_NT^2)
    # SHB = SHB/sum(SHB)
    return SHB
end
