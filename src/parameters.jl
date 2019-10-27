mutable struct Params{T<:Real}
    radius::T   # radius in cm
    pump_radius::T
    L::T        # length in cm
    L_eff::T    # effective length in cm
    h::T        # Planck constant
    c::T        # speed of light in m/s
    ev::T       # in J
    kB::T       # in eV/K
    T::T        # in K
    T_vA::AbstractVector
    T_vE::AbstractVector
    kBT::T      # in cm^-1
    M::T        # molecular mass in AMU
    norm_time::T# normalization to micro second

    σ_GKC::T    # in Angstrom^2, gas kinetic collision cross section
    σ_DD::T     # dipole - dipole
    σ_SPT::T
    σ_36::T     # v3->v6 cross section
    σ_VS::T     # V swap cross section

    v_avg::T    # in m/sec, average relative velocity between molecules
    kvs::T      # in m^3/microsec
    # kvsplit::T

    EG::T
    E3::T
    E6::T
    E23::T
    E36::T
    E26::T
    Q::T
    f_G_0::T
    f_3_0::T
    f_6_0::T
    f_GA::AbstractVector
    f_3A::AbstractVector
    f_6A::AbstractVector
    f_GE::AbstractVector
    f_3E::AbstractVector
    f_6E::AbstractVector

    f_23::T
    f_36::T
    f_26::T

    JL::Integer #pump: JL in V0 ->JU in V3
    JU::Integer
    K0::Integer

    pumpbranch::AbstractString #'R' or 'P' or 'Q'

    CL::T
    CL1::T
    CU::T
    CU1::T

    g_L::T
    g_U::T

    NA::T

    f₀::T
    f_offset::T
    f_pump::T
    f_dir_lasing::T
    f_ref_lasing::T
    Δ_f₀D::T
    f_range::T

    n_rot::Integer
    n_vib::Integer

    mu0::T
    eps0::T
    # mode_num: 1: TE01 / 2: TE12 / 3: TE02 / 4: TE22 / 5: TE11 / 6: TE21 / 7: TM01 / 8: TM11
    mode_num::Integer
    # zero point of the Bessel function:
    p_library::AbstractVector
    n0::T #refractive index
    t_spont::T
    Δν_THz::T

    ### parameters that are not system properties, that are related to pressure and power etc :
    pressure::T     # in unit mTorr
    power::T        # in unit W
    powerF::AbstractVector
    powerB::AbstractVector
    averagePF::AbstractVector
    averagePB::AbstractVector

    num_layers::Integer

    ntotal::T # in unit m^-3

    k63A::AbstractVector
    k36A::AbstractVector
    k63E::AbstractVector
    k36E::AbstractVector

    k3623::T
    k2336::T
    k2636::T
    k3626::T
    kro::T

    Δ_fP::T
    Δ_f_RabiF::AbstractVector
    Δ_f_RabiB::AbstractVector
    Δ_f_NTF::AbstractVector
    Δ_f_NTB::AbstractVector

    num_freq::Integer
    layer_unknown::Integer
    df::T
    f_dist_end::AbstractVector
    f_dist_ctr::AbstractVector
    f_dist_ctrB::AbstractVector

    velocity::AbstractVector
    f_dist_dir_lasing::AbstractVector
    f_dist_ref_lasing::AbstractVector

    gauss_dist::AbstractVector
    SHBF::AbstractArray
    SHBB::AbstractArray
    fp_lasing::AbstractVector
    fp_ref_lasing::AbstractVector

    alpha_0::T
    alpha_r::AbstractVector
    pump_0::T
    pump_IR::AbstractArray

    Δr::T
    r_int::AbstractVector

    kwall::AbstractVector
    MFP::T
    D::T

    kDDmat::AbstractArray
    ka::AbstractVector

    niter::Integer
    lin_solver::AbstractString
    model_flag::Integer
    solstart_flag::Integer

    backward::Integer
    ACStark::Integer
    effectiveL::Integer
    MultRef::Integer

    D_factor::T

    f_dirgain_dist::AbstractVector
    f_refgain_dist::AbstractVector
    dirgain::AbstractVector
    refgain::AbstractVector

    script::Integer

    evol_t::AbstractVector

    err_tv::Bool

    beta13::T

    # Wi::T
    WiU::T
    WiL::T

    optcavity::Bool

end

function Params(DefaultT=Float64;
    radius = 0.25, # in cm
    L = 15, # in cm
    T = 300,
    pressure = 100.0,
    power = 10.0,
    ####################################
    ## molecule setup
    ####################################
    M = 35,
    ####################################
    ## pump setup
    ####################################
    JL = 4,
    K0 = 3,
    pumpbranch = "R",
    f_offset = 30e6,
    n_rot = 18,
    n0 = 1.0,
    t_spont = 10,
    ####################################
    ## model/solver setup
    ####################################
    num_layers = 10,
    model_flag = 2,
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

    evol_t = 0:.1:1
    err_tv = false

    pump_radius = radius # deprecated
    L_eff = L

    ###################################################
    #### molecule setup
    ###################################################
    if M == 35 # 13CH3F
        σ_GKC = 44
        σ_DD = 320
        σ_SPT = 137
        σ_36 = 1.61
        σ_VS = 21
        EG = 0
        E3 = 1027.49
        E6 = 1175
        E23 = 2100
        E36 = 2250
        E26 = 2400
        A0, B0, DJ0, DJK0 = (155.3528, 24.8626427, 5.7683E-05, 0.00042441)
        A3, B3, DJ3, DJK3 = (155.3528, 24.5421324, 0.000055156, 0.00047788)
    elseif M == 34
        σ_GKC = 44
        σ_DD = 320
        σ_SPT = 50.0
        σ_36 = 3.21
        σ_VS = 18.9
        EG = 0
        E3 = 1048.6107008736
        E6 = 1182.35
        E23 = 2100
        E36 = 2250
        E26 = 2400
        A0, B0, DJ0, DJK0 = (155.3528, 25.5361499, 0.000060233, 0.0004395743) # in GHz
        A3, B3, DJ3, DJK3 = (155.3528, 25.1975092, 5.68788E-05, 0.000518083) # in GHz
    end
    g_6 = 2

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
    CL, CL1, CU, CU1 = compCs(JL, JU, K0, h, T, M)
    # println(E3, ", ", B3, ", ", DJ3, ", ", DJK3)
    EU = E3*c/1e7 + B3*JU*(JU+1) + (A3-B3)*K0^2 - DJ3*JU^2*(JU+1)^2 -DJK3*JU*(JU+1)*K0^2 # in GHz
    # println(EU)
    EL = B0*JL*(JL+1) + (A0-B0)*K0^2 - DJ0*JL^2*(JL+1)^2 -DJK0*JL*(JL+1)*K0^2
    # println(EL)
    f₀ = (EU - EL) * 1e9 # in Hz
    f_pump = f₀ + f_offset

    f_dir_lasing = ΔEr(JU-1, JU, K0, "V3", M)
    f_ref_lasing = ΔEr(JL, JL+1, K0, "V0", M)

    ###################################################
    #### derived molecular parameters
    ###################################################
    v_avg = 205*sqrt(T/M)
    vel = v_avg/sqrt(2)/norm_time # avg absolute vel in m/microsec
    kvs = v_avg*σ_VS * (1e-10)^2/norm_time
    MFP = 0.732*T/pressure/σ_GKC # in cm
    # σ_VSplit = 12.4
    # kvsplit = v_avg*σ_VSplit * (1e-10)^2/norm_time
    # diffusion coefficient in m^2/microsec. Einstein-Smoluchowski equation;
    D = 1/3 * vel * MFP * 1e-2 * D_factor
    kDD = 19.8 * pressure * σ_DD/ sqrt(T*M)

    Δ_f₀D = 3.58e-7*f₀*sqrt(T/M)
    Δ_fP = 15e6*(pressure/1e3)

    ntotal = 9.66e24 * pressure * 1e-3 / T # with unit m^-3
    kro = ntotal * v_avg * σ_VS * (1e-10)^2 / norm_time

    if model_flag==1
        Q = exp(EG) + exp(-E3/kBT) + 2*exp(-E6/kBT) + exp(-E23/kBT) +
            exp(-E36/kBT) + exp(-E26/kBT)
        f_G_0 = exp(EG)/Q
        f_3_0 = exp(-E3/kBT)/Q
        f_6_0 = g_6*exp(-E6/kBT)/Q
    elseif model_flag==2
        Q = Qv(kB, T)
        f_G_0 = exp(EG)/Q
        f_3_0 = exp(-E3/kBT)/Q
        f_6_0 = 1- f_G_0 - f_3_0
    end
    f_23 = exp(-E23/kBT)/Q
    f_36 = exp(-E36/kBT)/Q
    f_26 = exp(-E26/kBT)/Q

    if model_flag == 1
        k63 = ntotal*v_avg*σ_36*(1e-10)^2/norm_time/2 # in 1/microsec
        k36 = exp(-(E6-E3)/kBT) * k63 * g_6
    else
        k36 = 100.0
        k63 = k36 * f_3_0/f_6_0
    end
    k3623 = ntotal*v_avg*σ_36*(1e-10)^2/norm_time
    k2336 = exp(-(E36-E23)/kBT) * k3623
    k2636 = ntotal*v_avg*σ_36*(1e-10)^2/norm_time
    k3626 = exp(-(E26-E36)/kBT) * k2636

    beta13 = 1.20 * sqrt(power)/radius

    totalNL0 = CL * ntotal * f_G_0/2
    totalNU0 = CU * ntotal * f_3_0/2

    ###################################################
    #### solver/discretization parameters
    ###################################################
    Δr = radius/100 / num_layers # in m
    r_ext = linspace(0,radius/100,num_layers+1)
    r_int = 0.5*(r_ext[1:end-1] + r_ext[2:end]) # in m

    f_range = 5*Δ_f₀D
    num_freq = round(Int64,max(50,2*f_range/(Δ_fP/4)))
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

    # assign J levels
    if JU<=K0 || JL<K0 # JL>=K0-1+n_rot÷2 || JU>=K0+n_rot÷2 ||
        throw(ArgumentError("JL or JU is out of bounds!"))
    end

    if model_flag == 1
        n_vib = 12
    elseif model_flag == 2
        n_vib = 6
    end

    n_rot = 2 * (5+max(JU, JL+1)-K0+1)
    J = Jlevels(n_rot, JL, JU, K0) #3:(n_rot÷2+2)

    # in 1/microsec. rate_ij = rate_DD * prob. of ij collision, also kDDmat[i,j]:
    kDDmat = zeros(n_rot, n_rot)
    for i in 2:length(J)
        kDDmat[i, i-1] = kDDmat[i+n_rot÷2, i+n_rot÷2-1] =
        kDD*Q_selectn_hl(J[i], K0)/(1e3)
    end
    for i in 1:length(J)-1
        kDDmat[i, i+1] = kDDmat[i+n_rot÷2, i+n_rot÷2+1] =
        kDD*Q_selectn_lh(J[i], K0, M)/(1e3)
    end
    # K-swap rates -> goes to thermal pool, in 1/microsec
    ka = 19.8*pressure*σ_SPT/sqrt(T*M)/(1e3) * ones(n_rot)

    layer_unknown = n_rot*num_freq + n_vib

    ###################################################
    #### physical terms initialization
    ###################################################
    T_vA = T * ones(num_layers)
    T_vE = T * ones(num_layers)
    f_GA = f_G_0 * ones(num_layers)
    f_3A = f_3_0 * ones(num_layers)
    f_6A = f_6_0 * ones(num_layers)
    f_GE = f_G_0 * ones(num_layers)
    f_3E = f_3_0 * ones(num_layers)
    f_6E = f_6_0 * ones(num_layers)

    k63A = k63 * ones(num_layers)
    k63E = k63 * ones(num_layers)
    k36A = k36 * ones(num_layers)
    k36E = k36 * ones(num_layers)

    fp_lasing = f_NT_normalized(f_dist_dir_lasing, Δ_fP, f_dir_lasing, df*f_dir_lasing/f₀) # not used
    fp_ref_lasing = f_NT_normalized(f_dist_ref_lasing, Δ_fP, f_ref_lasing, df*f_ref_lasing/f₀) # not used

    f_dirgain_dist = linspace(f_dir_lasing-40e6, f_dir_lasing+40e6, 500)
    f_refgain_dist = linspace(f_ref_lasing-40e6, f_ref_lasing+40e6, 500)
    dirgain = zeros(size(f_dirgain_dist))
    refgain = zeros(size(f_refgain_dist))

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
    alpha_0 = exp(-log(2)*((f_pump-f₀)/Δ_f₀D)^2)*sqrt(log(2)/pi)/Δ_f₀D *
              8*pi^3/3/h/c * (totalNL0 - totalNU0) *
              1e-36 * (0.2756^2*16.0/45) * f_pump * 1e-13
    alpha_r = alpha_0 * ones(num_layers)
    pump_IR = zeros(num_freq, num_layers)
    # in unit m^-3 microsec^-1
    pump_0 = 9.4e13 * power/(radius^2)/Δ_f₀D * (0.2756^2*16.0/45) *
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


function lorentz_dist(ν, Δ_f_NT, f_pump)
    return  1/π .* Δ_f_NT ./ ((ν - f_pump).^2 + Δ_f_NT^2)
end

function f_NT_ampl(ν, Δ_f_NT, f_pump)
    SHB = Δ_f_NT^2 ./ ((ν - f_pump).^2 + Δ_f_NT^2)
    # SHB = SHB/sum(SHB)
    return SHB
end

function f_NT_normalized(ν, Δ_f_NT, f_pump, df)
    SHB = 1/π * df * Δ_f_NT ./ ((ν .- f_pump).^2 .+ Δ_f_NT^2)
    # SHB = SHB/sum(SHB)
    return SHB
end

function emission_broaden(ν, vi, p, df)
  # emission spectrum of subclass vi at frequencys nu
  spectrum = zeros(size(ν))
  f31 = p.f_dist_ctr[vi]
  f32 = p.f_dist_dir_lasing[vi]
  deltap = 2π * (p.f_pump - f31)
  # deltap = 2π * 25e6
  γ = 0.5 * sqrt(deltap^2 + 4*p.beta13^2)
  τ = 1/(2π*p.Δ_fP)
  for i in 1:length(ν)
    Ω = 0.5*deltap - 2π * (ν[i]-f32)
    spectrum[i] = 1/(1+(γ-Ω)^2*τ^2) + 1/(1+(γ+Ω)^2*τ^2) +
    (2*(γ^2-Ω^2)*τ^2*(1+2*γ^2*τ^2)-2)/(1+4γ^2*τ^2)/(1+(γ-Ω)^2*τ^2)/(1+(γ+Ω)^2*τ^2)
  end
  spectrum = df*(1+4*γ^2*τ^2)/(4*γ^2*τ) * spectrum
  return spectrum
end



function ΔEr(J1, J2, K, V, M) # in Hz
    if V=="V0" && M==35
        B = 24862.6427e6
        DJ = 0.057683e6
        DJK = 0.42441e6
    elseif V=="V3" && M==35
        B = 24542.1324e6
        DJ = 0.055156e6
        DJK = 0.47788e6
    elseif V=="V0" && M==34
        B = 25.5361499e9
        DJ = 0.000060233e9
        DJK = 0.000439574e9
    elseif V=="V3" && M==34
        B = 25.1975092e9
        DJ = 5.68788E4
        DJK = 0.000518083e9
    end
    return (B-DJK*K^2)*(J2*(J2+1)-J1*(J1+1)) - DJ*(-J1^2*(J1+1)^2+J2^2*(J2+1)^2)
end

function derv_bessel(ν,x)
    return 0.5.*(besselj(ν-1,x)-besselj(ν+1,x))
end

function compCs(JL, JU, K0, h, T, M)
    JKE = JKEgenerate(100, M) # K from -J to J
    Qi, Q, QA = Qcompute(JKE, h, T)
    CL = CL1 = CU = CU1 = 0.
    for i in 1:size(JKE, 1)
        if JKE[i, 1] == JL && JKE[i, 2] == K0
            CL = K0==0 ? Qi[i]/QA : Qi[i]/QA * 2 # * 2 because ±K0
        elseif JKE[i, 1] == JL+1 && JKE[i, 2] == K0
            CL1 = K0==0 ? Qi[i]/QA : Qi[i]/QA * 2
        end
    end

    for i in 1:size(JKE, 1)
        if JKE[i, 1] == JU && JKE[i, 2] == K0
            CU = K0==0 ? Qi[i]/QA : Qi[i]/QA * 2
        elseif JKE[i, 1] == JU-1 && JKE[i, 2] == K0
            CU1 = K0==0 ? Qi[i]/QA : Qi[i]/QA * 2
        end
    end

    return (CL, CL1, CU, CU1)
end

function Qcompute(JKE, h, T)
    # Q: total partition func
    # QA: partition func of A type, same with E type
    # Qi: partition func of (J, K), K can be negative
    J = JKE[:,1]
    K = JKE[:,2]
    E = JKE[:,3]
    n = length(J)
    Q = QA = 0.0
    Qi = zeros(n)
    for i = 1:n
        gi = 2J[i] + 1
        if K[i]%3 == 0 # A type
            gi *= 2
        end
        qi = gi * exp(-E[i]*1e9*h/(1.38064852e-23*T))
        Q += qi
        if K[i]%3 == 0
            QA += qi
        end
        Qi[i] = qi
    end
    return (Qi, Q, QA)
end

function JKEgenerate(n, M)
    # generate J, K, E array up to J = n
    if M == 35 # 13CH3F
        A, B, DJ, DJK = (155.3528, 24.8626427, 5.7683E-05, 0.00042441)
    elseif M == 34 # 12CH3F
        A, B, DJ, DJK = (155.3528, 25.5361499, 0.000060233, 0.000439574)
    else
        throw(ArgumentError("M can only be 34 (12CH3F) or 35 (13CH3F)"))
    end
    JKE = zeros((n+1)^2, 3)
    linenum = 0
    for j in 0:n
        for k in -j:1:j
            linenum += 1
            JKE[linenum, 1] = j
            JKE[linenum, 2] = k
            JKE[linenum, 3] = B*j*(j+1) + (A-B)*k^2 - DJ*j^2*(j+1)^2 -DJK*j*(j+1)*k^2
        end
    end
    return JKE
end

function Jlevels(n_rot, JL, JU, K0)
    Jmin = K0 # max(K0, min(JL-n_rot÷4, JU-n_rot÷4))
    Jmax = n_rot÷2 + Jmin - 1
    if Jmax < JU+2 || Jmax < JL+2
        throw(ArgumentError("number of rotational levels need to be increased!"))
    end
    return Jmin:Jmax
end
