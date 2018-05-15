
function CO(DefaultT=Float64;
    radius = 0.25, # in cm
    L = 15, # in cm
    T = 300,
    pressure = 100.0,
    power = 10.0,
    ####################################
    ## molecule setup
    ####################################
    M = 28,
    ####################################
    ## pump setup
    ####################################
    JL = 4,
    pumpbranch = "R",
    f_offset = 30e6,
    n_rot = 18,
    n0 = 1.0,
    t_spont = 10,
    ####################################
    ## model/solver setup
    ####################################
    num_layers = 10,
    solstart_flag = 0,
    optcavity = false,
    D_factor = 1.0,
    WiU = 0,
    WiL = 0
    )

    ###################################################
    #### physical constants
    ###################################################
    h = 6.626068e-34
    c = 2.99792458e8
    ev = 1.60217646e-19
    kB = 8.617342e-5 # in ev/K
    mu0 = 4e-7*pi
    eps0 = 8.85e-12
    NA = 6.0221413e23
    kBT = kB*T*8065.73 # in cm^-1

    #### unit constants
    norm_time = 1e6

    ###################################################
    #### model configuration
    ###################################################
    niter = 10 # deprecated
    lin_solver = "Default" # deprecated
    backward = 1 # not used
    ACStark = 1 # not used
    effectiveL = 0 # deprecated
    script = 0 # deprecated
    MultRef = 1 # deprecated
    Δν_THz = 25e6 # deprecated
    mode_num = 1 # deprecated
    # mode_num: 1: TE01 / 2: TE12 / 3: TE02 / 4: TE22 / 5: TE11 / 6: TE21 / 7: TM01 / 8: TM11
    #zeros of Bessel functions:
    p_library = [3.83, 5.33, 7.02, 6.71, 1.84, 3.05, 2.4, 3.83] # deprecated
    K0 = 0
    model_flag = 1
    σ_SPT = 0
    σ_36 = 0
    σ_VS = 0
    E6 = 4260
    E23 = 6350
    E36 = 8415
    E26 = 0


    evol_t = 0:.1:1
    err_tv = false

    pump_radius = radius # deprecated
    L_eff = L

    ###################################################
    #### molecule setup
    ###################################################
    σ_GKC = 44
    σ_DD = 3.15*sqrt(T*M)/3.15 * 3.2
    EG = 0
    E3 = 2143.0
    (B, DJ) = (57.635968, 0.1835058E-6) # in GHz

    ###################################################
    #### pump setup
    ###################################################

    if pumpbranch == "R"
        JU = JL+1
    elseif pumpbranch == "P"
        JU = JL-1
    else
        throw(ArgumentError("pump branch can only be P or R!"))
    end
    g_L = 2JL + 1
    g_U = 2JU + 1
    CL, CL1, CU, CU1 = compCsCO(JL, JU, h, T, M)

    EU = E3*c/1e7 + B*JU*(JU+1) - DJ*JU^2*(JU+1)^2 # in GHz
    EL = B*JL*(JL+1) - DJ*JL^2*(JL+1)^2
    f₀ = (EU - EL) * 1e9 # in Hz
    f_pump = f₀ + f_offset

    f_dir_lasing = 2B*JU - 4DJ*JU^3
    f_ref_lasing = 2B*(JL+1) - 4DJ*(JL+1)^3

    ###################################################
    #### derived molecular parameters
    ###################################################
    v_avg = 205*sqrt(T/M)
    vel = v_avg/sqrt(2)/norm_time # avg absolute vel in m/microsec
    kvs = 0
    MFP = 0.732*T/pressure/σ_GKC # in cm
    # diffusion coefficient in m^2/microsec. Einstein-Smoluchowski equation;
    D = 1/3 * vel * MFP * 1e-2 * D_factor
    kDD = 19.8 * pressure * σ_DD/ sqrt(T*M)

    Δ_f₀D = 3.58e-7*f₀*sqrt(T/M)
    Δ_fP = 3.2e6*(pressure/1e3)

    ntotal = 9.66e24 * pressure * 1e-3 / T # with unit m^-3
    kro = 0

    Q = exp(EG) + exp(-E3/kBT)
    f_G_0 = exp(EG)/Q
    f_3_0 = exp(-E3/kBT)/Q

    f_6_0 = f_23 = f_36 = f_26 = 0.0

    k63 = k36 = 0.0
    k3623 = k2336 = k2636 = k3626 = 0.0

    beta13 = 1.20 * sqrt(power)/radius

    ###################################################
    #### solver/discretization parameters
    ###################################################
    Δr = radius/100 / num_layers # in m
    r_ext = linspace(0,radius/100,num_layers+1)
    r_int = 0.5*(r_ext[1:end-1] + r_ext[2:end]) # in m

    f_range = 2*Δ_f₀D
    num_freq = round(Int64,max(50,f_range/(Δ_fP/4)))
    num_freq = 1
    df = 2.0 * f_range / num_freq
    f_dist_end = linspace(-f_range, f_range, num_freq + 1) + f₀
    f_dist_ctr = f_dist_end[1:end-1] + df/2

    velocity = (f_dist_ctr - f₀)/f₀ # in unit c
    f_dist_dir_lasing = velocity * f_dir_lasing + f_dir_lasing
    f_dist_ref_lasing = velocity * f_ref_lasing + f_ref_lasing

    f_dist_ctrB = f₀ - f₀ * velocity

    norm_dist = Normal(f₀, Δ_f₀D / sqrt(2*log(2)))
    pdf1 = pdf(norm_dist, f_dist_ctr)
    gauss_dist = pdf1 / sum(pdf1)
    # gauss_dist = pdf1 * df

    ## n_rot
    n_vib = 0

    # n_rot = 2 * (5+max(JU, JL+1))
    J = 0:(n_rot÷2-1)

    # in 1/microsec. rate_ij = rate_DD * prob. of ij collision, also kDDmat[i,j]:
    kDDmat = zeros(n_rot, n_rot)
    for i in 2:length(J)
        kDDmat[i, i-1] = kDDmat[i+n_rot÷2, i+n_rot÷2-1] =
        kDD*Q_selectn_hl(J[i], 0)/(1e3)
    end
    for i in 1:length(J)-1
        kDDmat[i, i+1] = kDDmat[i+n_rot÷2, i+n_rot÷2+1] =
        kDD*Q_selectn_lh(J[i], 0)/(1e3)
    end
    # K-swap rates -> goes to thermal pool, in 1/microsec
    ka = [0.0]

    layer_unknown = n_rot*num_freq

    ###################################################
    #### physical terms initialization
    ###################################################
    T_vA = T_vE = [0]
    f_GA = f_3A = f_6A = f_GE = f_3E = f_6E = [0]
    k63A = k63E = k36A = k36E = [0]

    fp_lasing = f_NT_normalized(f_dist_dir_lasing, Δ_fP, f_dir_lasing, df*f_dir_lasing/f₀) # not used
    fp_ref_lasing = f_NT_normalized(f_dist_ref_lasing, Δ_fP, f_ref_lasing, df*f_ref_lasing/f₀) # not used

    f_dirgain_dist = linspace(f_dir_lasing-40e6, f_dir_lasing+40e6, 500) # not used
    f_refgain_dist = linspace(f_ref_lasing-40e6, f_ref_lasing+40e6, 500) # not used
    dirgain = zeros(size(f_dirgain_dist)) # not used
    refgain = zeros(size(f_refgain_dist)) # not used

    Δ_f_RabiF = zeros(num_layers)
    Δ_f_RabiB = zeros(num_layers)
    Δ_f_NTF = ones(num_layers) * Δ_fP
    Δ_f_NTB = ones(num_layers) * Δ_fP
    SHBF = zeros(num_freq, num_layers)
    SHBB = zeros(num_freq, num_layers)
    powerF = zeros(num_layers)
    powerB = zeros(num_layers)
    averagePF = power * ones(num_layers)
    averagePB = power * ones(num_layers)

    # alpha_0 in m^-1
    totalNL0 = CL * ntotal * f_G_0
    totalNU0 = CU * ntotal * f_3_0
    alpha_0 = exp(-log(2)*((f_pump-f₀)/Δ_f₀D)^2)*sqrt(log(2)/pi)/Δ_f₀D *
              8*pi^3/3/h/c * (totalNL0 - totalNU0) *
              1e-36 * (0.01156806*max(JL, JU)/(2JL+1)) * f_pump * 1e-13
    alpha_r = alpha_0 * ones(num_layers)
    pump_IR = zeros(num_freq, num_layers)
    # in unit m^-3 microsec^-1
    pump_0 = 9.4e13 * power/(radius^2)/Δ_f₀D * (0.01156806*max(JL, JU)/(2JL+1)) *
                exp(-log(2)*((f_pump-f₀)/Δ_f₀D)^2)/norm_time

    # kwall = WallRate(radius, pressure, r_int, ntotal, M, T, NA, v_avg, σ_GKC) + 1e-10
    kwall = zeros(num_layers)

    return Params{DefaultT}(radius, pump_radius, L, L_eff, h, c, ev, kB, T, T_vA, T_vE, kBT, M, norm_time,
    σ_GKC, σ_DD, σ_SPT, σ_36, σ_VS,
    v_avg, kvs, #kvsplit,
    EG, E3, E6, E23, E36, E26, Q,
    f_G_0, f_3_0, f_6_0, f_GA, f_3A, f_6A, f_GE, f_3E, f_6E, f_23, f_36, f_26,
    JL, JU, K0, pumpbranch, CL, CL1, CU, CU1, g_L, g_U, NA,
    f₀, f_offset, f_pump, f_dir_lasing, f_ref_lasing, Δ_f₀D, f_range,
    n_rot, n_vib,
    mu0, eps0,
    mode_num, p_library, n0, t_spont, Δν_THz,
    pressure, power, powerF, powerB, averagePF, averagePB,
    num_layers, ntotal,
    k63A, k36A, k63E, k36E, k3623, k2336, k2636, k3626, kro,
    Δ_fP, Δ_f_RabiF, Δ_f_RabiB, Δ_f_NTF, Δ_f_NTB,
    num_freq, layer_unknown, df, f_dist_end, f_dist_ctr, f_dist_ctrB,
    velocity, f_dist_dir_lasing, f_dist_ref_lasing,
    gauss_dist, SHBF, SHBB, fp_lasing, fp_ref_lasing,
    alpha_0, alpha_r, pump_0, pump_IR,
    Δr, r_int,
    kwall, MFP, D,
    kDDmat,
    ka,
    niter, lin_solver, model_flag, solstart_flag,
    backward, ACStark, effectiveL, MultRef, D_factor,
    f_dirgain_dist, f_refgain_dist, dirgain, refgain,
    script,
    evol_t,
    err_tv,
    beta13,
    WiU, WiL,
    optcavity
    )
end

function compCsCO(JL, JU, h, T, M)
    cL = compCCO(JL, h, T, M)
    cL1 = compCCO(JL+1, h, T, M)
    cU = compCCO(JU, h, T, M)
    cU1 = compCCO(JU-1, h, T, M)
    return (cL, cL1, cU, cU1)
end

function compCCO(JL, h, T, M)
    # compute fraction of JL in that vibrational level
    Js = 0:1:100 # K from -J to J
    Q = ql = 0.0
    (B, DJ) = (57.635968, 0.1835058E-6) # in GHz
    for J in Js
        Ej = B*J*(J+1) - DJ*J^2*(J+1)^2
        qj = (2J+1) * exp(-Ej*1e9*h/(1.38064852e-23*T))
        Q += qj
        if J == JL
            ql = qj
        end
    end
    return ql/Q
end
