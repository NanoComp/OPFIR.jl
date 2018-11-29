function N2O(DefaultT=Float64;
    radius = 0.25, # in cm
    L = 15, # in cm
    T = 300,
    pressure = 100.0,
    power = 10.0,
    ####################################
    ## molecule setup
    ####################################
    M = 44,
    ####################################
    ## pump setup
    ####################################
    JL = 4,
    pumpbranch = "R",
    f_offset = 0.e6,
    n_rot = 18,
    n0 = 1.0,
    t_spont = 10,
    ####################################
    ## model/solver setup
    ####################################
    num_layers = 8,
    solstart_flag = 0,
    optcavity = false,
    D_factor = 1.0,
    WiU = 0,
    WiL = 0,
    # approach = 1: approach 1 only SPT, between rot levels;
    # approach = 2: approach 2 SPT, VV, direct from vib levels;
    # approach = 3: vib thermal pools modeling, suggested by Henry;
    approach = 1,
    n_vib = 1
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
    backward = 1 # not used, assumed to be true
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
    model_flag = 2 # default to be 2 for N2O; 1 for CO;
    σ_SPT = 0
    σ_36 = 0
    σ_VS = 0
    E6 = 0
    E23 = 0
    E36 = 0
    E26 = 0

    evol_t = 0:.1:1
    err_tv = false

    pump_radius = radius # deprecated
    L_eff = L

    ###################################################
    #### molecule setup
    ###################################################
    σ_GKC = 15
    # σ_DD = sqrt(T*M)/3.15 * 4.0 # ~146 A^2
    σ_DD = 146

    EG = 0
    E3 = 2224.0
    B, DJ = [0.419, 17.6e-8]*c/1e7 # in GHz

    ###################################################
    #### pump setup
    ###################################################

    if pumpbranch == "R"
        JU = JL+1
    elseif pumpbranch == "Q"
        JU = JL
    elseif pumpbranch == "P"
        JU = JL-1
    else
        throw(ArgumentError("pump branch can only be P, Q, R!"))
    end
    g_L = 2JL + 1
    g_U = 2JU + 1
    CL, CL1, CU, CU1 = compCsN2O(JL, JU, h, T, M)
    # println(E3, ", ", B3, ", ", DJ3, ", ", DJK3)
    EU = E3*c/1e7 + B*JU*(JU+1) - DJ*JU^2*(JU+1)^2 # in GHz
    EL = B*JL*(JL+1) - DJ*JL^2*(JL+1)^2
    f₀ = (EU - EL) * 1e9 # in Hz
    f_pump = f₀ + f_offset

    f_dir_lasing = (2B*JU - 4DJ*JU^3) * 1e9
    f_ref_lasing = (2B*(JL+1) - 4DJ*(JL+1)^3) * 1e9

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
    Δ_fP = 4.0e6*(pressure/1e3)

    ntotal = 9.66e24 * pressure * 1e-3 / T # with unit m^-3
    kro = 0

    Q = Qv(kB, T)
    f_G_0 = exp(EG)/Q
    f_3_0 = exp(-E3/kBT)/Q
    f_6_0 = 1- f_G_0 - f_3_0

    f_23 = f_36 = f_26 = 0.0

    k63 = k36 = 0.0
    k3623 = k2336 = k2636 = k3626 = 0.0

    beta13 = 0 * sqrt(power)/radius

    ###################################################
    #### solver/discretization parameters
    ###################################################
    Δr = radius/100 / num_layers # in m
    r_ext = linspace(0,radius/100,num_layers+1)
    r_int = 0.5*(r_ext[1:end-1] + r_ext[2:end]) # in m

    # f_range = 2*Δ_f₀D
    f_range = 80Δ_fP
    num_freq = round(Int64,max(50,2f_range/(Δ_fP/4)))

    df = 2.0 * f_range / num_freq
    f_dist_end = linspace(-f_range, f_range, num_freq + 1) + f₀
    f_dist_ctr = f_dist_end[1:end-1] + df/2

    velocity = (f_dist_ctr - f₀)/f₀ # in unit c
    f_dist_dir_lasing = velocity * f_dir_lasing + f_dir_lasing
    f_dist_ref_lasing = velocity * f_ref_lasing + f_ref_lasing

    f_dist_ctrB = f₀ - f₀ * velocity

    norm_dist = Normal(f₀, Δ_f₀D / sqrt(2*log(2)))
    pdf1 = pdf.(norm_dist, f_dist_ctr)
    gauss_dist = pdf1 / sum(pdf1)

    J, n_rot = JlevelsN2O(JL, JU, K0)

    # in 1/microsec. rate_ij = rate_DD * prob. of ij collision, also kDDmat[i,j]:
    kDDmat = zeros(n_rot, n_rot)
    for i in 2:length(J)
        kDDmat[i, i-1] = kDDmat[i+n_rot÷2, i+n_rot÷2-1] =
        kDD*Q_selectn_hl(J[i], 0)/(1e3)
    end
    for i in 1:length(J)-1
        dE = (2B*(J[i]+1) - 4DJ*(J[i]+1)^3) * 1e9 # in Hz
        kDDmat[i, i+1] = kDDmat[i+n_rot÷2, i+n_rot÷2+1] =
        kDD*Q_selectn_lh(J[i], K0, M)/(1e3) * exp(-dE * h/1.38e-23/T)
    end

    k_rottherm0 = 19.8 * pressure * σ_GKC/ sqrt(T*M) # in msec-1
    if approach == 2
        ks = kGKC(n_rot, J, h, kB, T, M) * k_rottherm0/1000 # in microsec-1
        kDDmat += ks
    elseif approach == 1
        k_rottherm = zeros(n_rot÷2, n_rot÷2)
        for si in 1:length(J)
            statei = J[si]
            for sf in 1:length(J)
                statef = J[sf]
                cf = compCN2O(statef, h, T, M)
                k_rottherm[si, sf] = k_rottherm0 * cf /1000 # in 1/microsec

                kDDmat[si, sf] += k_rottherm[si, sf]
                kDDmat[si+n_rot÷2, sf+n_rot÷2] += k_rottherm[si, sf]
            end
        end
    end
    
    # K-SPT
    ka = [k_rottherm0/1000.0] # in microsec
    ## rotational population fraction to all
    rotpopfr = rotpopfracl(h, T, M, n_rot, f_G_0, f_3_0, J)
    ck = zeros(n_rot÷2)
    for k in 1:n_rot÷2
        ck[k] = compCN2O(J[k], h, T, M)
    end
    ck = ck/sum(ck)
    cj = vcat(ck, ck)

    layer_unknown = n_rot*num_freq + n_vib

    ###################################################
    #### physical terms initialization
    ###################################################
    T_vA = T * ones(num_layers)
    T_vE = T * ones(num_layers) # deprecated for N2O
    f_GA = f_G_0 * ones(num_layers)
    f_3A = f_3_0 * ones(num_layers)
    f_6A = f_6_0 * ones(num_layers)
    f_GE = f_G_0 * ones(num_layers) # deprecated for N2O
    f_3E = f_3_0 * ones(num_layers) # deprecated for N2O
    f_6E = f_6_0 * ones(num_layers) # deprecated for N2O

    # ks are all 0
    k63A = k63 * ones(num_layers)
    k63E = k63 * ones(num_layers)
    k36A = k36 * ones(num_layers)
    k36E = k36 * ones(num_layers)

    ## not used anymore:
    fp_lasing = f_NT_normalized(f_dist_dir_lasing, Δ_fP, f_dir_lasing, df*f_dir_lasing/f₀) # not used
    fp_ref_lasing = f_NT_normalized(f_dist_ref_lasing, Δ_fP, f_ref_lasing, df*f_ref_lasing/f₀) # not used

    f_dirgain_dist = linspace(f_dir_lasing-40e6, f_dir_lasing+40e6, 500)
    f_refgain_dist = linspace(f_ref_lasing-40e6, f_ref_lasing+40e6, 500)
    dirgain = zeros(size(f_dirgain_dist))
    refgain = zeros(size(f_refgain_dist))

    ##
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
              1e-36 * (0.2457^2*max(JL, JU)/(2JL+1)) * f_pump * 1e-13

    alpha_r = zeros(num_layers)
    pump_IR = zeros(num_freq, num_layers)
    # in unit m^-3 microsec^-1
    pump_0 = 9.4e13 * power/(radius^2)/Δ_f₀D * (0.2457^2*max(JL, JU)/(2JL+1)) *
                exp(-log(2)*((f_pump-f₀)/Δ_f₀D)^2)/norm_time

    # kwall = WallRate(radius, pressure, r_int, ntotal, M, T, NA, v_avg, σ_GKC) + 1e-10
    kwall = zeros(num_layers)

    return ParamsN2O{DefaultT}(radius, pump_radius, L, L_eff, h, c, ev, kB, T, T_vA, T_vE, kBT, M, norm_time,
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
    optcavity,
    rotpopfr, cj,
    approach
    )
end


function compCsN2O(JL, JU, h, T, M)
    cL = compCN2O(JL, h, T, M)
    cL1 = compCN2O(JL+1, h, T, M)
    cU = compCN2O(JU, h, T, M)
    cU1 = compCN2O(JU-1, h, T, M)
    return (cL, cL1, cU, cU1)
end

function compCN2O(JL, h, T, M)
    # compute fraction of JL in that vibrational level
    Js = 0:1:100
    Q = ql = 0.0
    B, DJ = [0.419, 17.6e-8]*2.99792458e8/1e7 # in GHz
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



function Q_selectn_hl(J, K)
    return (J^2-K^2)/(J*(2*J+1))
end

function Q_selectn_lh(J, K, M)
    return ((J+1)^2-K^2)/((J+1)*(2*J+1))
end

function Qv(kB, T)
    data = viblevelsN2O()
    Q = 1.0
    for i in 1:size(data, 1)
        Q += data[i,2] * exp(-data[i, 1]/(kB*T*8065.73))
    end
    return Q
end

function rotpopfracl(h, T, M, n_rot, f_G_0, f_3_0, J)
    ctot = 0.0
    cj = zeros(n_rot)
    Js = vcat(J, J)
    for k in 1:n_rot
        cj[k] = compCN2O(Js[k], h, T, M)
        if k <= n_rot ÷ 2
            cj[k] *= f_G_0
        else
            cj[k] *= f_3_0
        end
        ctot += cj[k]
    end
    return cj/ctot * (f_G_0+f_3_0)
end


function JlevelsN2O(JL, JU, K0)
    N = max(16, JL)
    Jmin = 0 #max(K0, JL-N)
    Jmax = Jmin + 2N
    # if Jmax < JU+2 || Jmax < JL+2
    #     throw(ArgumentError("number of rotational levels need to be increased!"))
    # end
    return Jmin:Jmax, (Jmax-Jmin+1)*2
end


function krottherm(kB, T, Vi)
    data = viblevelsN2O() # doesn't include ground state 0
    if Vi == "V3"
        E0 = data[6,1]
    elseif Vi == "V0"
        E0 = 0.0
    end
    Q = exp(-abs(E0-0)/(kB*T*8065.73))
    for i in 1:size(data, 1)
        Q += data[i,2] * exp(-abs(E0-data[i, 1])/(kB*T*8065.73))
    end

    if Vi == "V3"
        c1 = exp(-abs(data[6,1])/(kB*T*8065.73))/Q
        c2 = 1/Q
        c3 = 1-c1-c2
    elseif Vi == "V0"
        c1 = 1/Q
        c2 = exp(-abs(data[6,1])/(kB*T*8065.73))/Q
        c3 = 1-c1-c2
    end
    return c1, c2, c3
end


function kGKC(n_rot, J, h, kB, T, M)
    ks = zeros(n_rot, n_rot)
    # V0 -> V0 and V3
    c00, c03, c0Σ = krottherm(kB, T, "V0")
    for i in 1:n_rot÷2
        for j in 1:n_rot÷2
            ks[i, j] = compCN2O(J[j], h, T, M) * c00
        end
        for j in n_rot÷2+1:n_rot
            ks[i, j] = compCN2O(J[j-n_rot÷2], h, T, M) *c03
        end
    end
    # V3 -> V0 and V3
    c30, c33, c3Σ = krottherm(kB, T, "V3")
    for i in n_rot÷2+1:n_rot
        for j in 1:n_rot÷2
            ks[i, j] = compCN2O(J[j], h, T, M) * c30
        end
        for j in n_rot÷2+1:n_rot
            ks[i, j] = compCN2O(J[j-n_rot÷2], h, T, M) * c33
        end
    end
    return ks
end