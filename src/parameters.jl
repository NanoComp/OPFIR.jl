include("wallrates.jl")

type Params{T<:Real}
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

    C3L::T
    C4L::T
    C5L::T
    C4U::T
    C5U::T

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
    n0::T
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

    netrate_36A::AbstractVector

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

    k98_G::T
    k87_G::T
    k76_G::T
    k65_G::T
    k54_G::T
    k43_G::T
    k32_G::T
    k21_G::T

    k89_G::T
    k78_G::T
    k67_G::T
    k56_G::T
    k45_G::T
    k34_G::T
    k23_G::T
    k12_G::T

    k98_3::T
    k87_3::T
    k76_3::T
    k65_3::T
    k54_3::T
    k43_3::T
    k32_3::T
    k21_3::T

    k89_3::T
    k78_3::T
    k67_3::T
    k56_3::T
    k45_3::T
    k34_3::T
    k23_3::T
    k12_3::T

    k1a::T
    k2a::T
    k3a::T
    k4a::T
    k5a::T
    k6a::T
    k7a::T
    k8a::T
    k9a::T
    k10a::T
    k11a::T
    k12a::T
    k13a::T
    k14a::T
    k15a::T
    k16a::T
    k17a::T
    k18a::T

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
    radius = 0.25,
    pump_radius = radius,
    L = 15,
    L_eff = 15,
    h = 6.626068e-34,
    c = 3.0e8,
    ev = 1.60217646e-19,
    kB = 8.617342e-5, # in ev/K
    T = 300,
    M = 35,
    norm_time = 1e6,
    σ_GKC = 44,
    σ_DD = 320,
    σ_SPT = 137,
    σ_36 = 1.61,
    σ_VS = 21,
    EG = 0,
    E3 = 1050,
    E6 = 1200,
    E23 = 2100,
    E36 = 2250,
    E26 = 2400,
    C3L = 0.005926302*2,
    C4L = 0.007374597*2,
    C5L = 0.008652703*2,
    C4U = C4L,
    C5U = C5L,
    g_L = 9.0,
    g_U = 11.0,
    NA = 6.0221413e23,
    f₀ = 31.042748176e12,
    f_offset = 30e6,
    f_pump = f₀ + f_offset,
    f_dir_lasing = 245.38e9,
    f_ref_lasing = 248.56e9,
    n_rot = 18,
    mu0 = 4e-7*pi,
    eps0 = 8.85e-12,
    # mode_num: 1: TE01 / 2: TE12 / 3: TE02 / 4: TE22 / 5: TE11 / 6: TE21 / 7: TM01 / 8: TM11
    mode_num = 1,
    #zeros of Bessel functions:
    p_library = [3.83, 5.33, 7.02, 6.71, 1.84, 3.05, 2.4, 3.83],
    n0 = 1.0,
    t_spont = 100,
    Δν_THz = 25e6,
    ## pressure and power related parameters:
    pressure = 100.0,
    power = 10.0,
    num_layers = 10,
    alpha_0 = 65.0,
    niter = 10,
    lin_solver = "Default",
    model_flag = 1,
    solstart_flag = 0,
    backward = 1,
    ACStark = 1,
    effectiveL = 0,
    MultRef = 1,
    D_factor = 1.0,
    script = 0,
    evol_t = 0:.1:1,
    err_tv = false,
    WiU = 0,
    WiL = 0,
    optcavity = false,
    )

    if model_flag == 1
        n_vib = 12
    elseif model_flag == 2
        n_vib = 6
    end

    T_vA = T * ones(num_layers)
    T_vE = T * ones(num_layers)

    kBT = kB*T*8065.73 # in cm^-1
    v_avg = 205*sqrt(T/M)
    kvs = v_avg*σ_VS * (1e-10)^2/norm_time

    # σ_VSplit = 12.4
    # kvsplit = v_avg*σ_VSplit * (1e-10)^2/norm_time

    g_6 = 2
    netrate_36A = zeros(num_layers)
    if model_flag==1
        Q = exp(EG) + exp(-E3/kBT) + 2*exp(-E6/kBT) + exp(-E23/kBT) +
            exp(-E36/kBT) + exp(-E26/kBT)
        f_G_0 = exp(EG)/Q
        f_3_0 = exp(-E3/kBT)/Q
        f_6_0 = g_6*exp(-E6/kBT)/Q
        f_GA = f_G_0 * ones(num_layers)
        f_3A = f_3_0 * ones(num_layers)
        f_6A = f_6_0 * ones(num_layers)
        f_GE = f_G_0 * ones(num_layers)
        f_3E = f_3_0 * ones(num_layers)
        f_6E = f_6_0 * ones(num_layers)
    elseif model_flag==2
        Q = Qv(kB, T, script)
        f_G_0 = exp(EG)/Q
        f_3_0 = exp(-E3/kBT)/Q
        f_6_0 = 1- f_G_0 - f_3_0
        f_GA = f_G_0 * ones(num_layers)
        f_3A = f_3_0 * ones(num_layers)
        f_6A = f_6_0 * ones(num_layers)
        f_GE = f_G_0 * ones(num_layers)
        f_3E = f_3_0 * ones(num_layers)
        f_6E = f_6_0 * ones(num_layers)
        # println(f_G_0, f_3_0, f_6_0)
    end
    f_23 = exp(-E23/kBT)/Q
    f_36 = exp(-E36/kBT)/Q
    f_26 = exp(-E26/kBT)/Q

    Δ_f₀D = 3.58e-7*f₀*sqrt(T/M)
    f_range = 5*Δ_f₀D

    ### parameters that are directly related to pressure and power:
    ntotal = 9.66e24 * pressure * 1e-3 / T # with unit m^-3

    # k63 = ntotal*v_avg*σ_36*(1e-10)^2/norm_time/2 # in 1/microsec
    # k36 = exp(-(E6-E3)/kBT) * k63 * g_6
    # k36 = f_6_0/f_3_0 * k63
    if model_flag == 1
        k63 = ntotal*v_avg*σ_36*(1e-10)^2/norm_time/2 # in 1/microsec
        k36 = exp(-(E6-E3)/kBT) * k63 * g_6
    else
        k36 = 100.0
        k63 = k36 * f_3_0/f_6_0
    end

    k63A = k63 * ones(num_layers)
    k63E = k63 * ones(num_layers)
    k36A = k36 * ones(num_layers)
    k36E = k36 * ones(num_layers)

    k3623 = ntotal*v_avg*σ_36*(1e-10)^2/norm_time
    k2336 = exp(-(E36-E23)/kBT) * k3623
    k2636 = ntotal*v_avg*σ_36*(1e-10)^2/norm_time
    k3626 = exp(-(E26-E36)/kBT) * k2636
    kro = ntotal * v_avg * σ_VS * (1e-10)^2 / norm_time
    # kro = 0

    Δ_fP = 15e6*(pressure/1e3)

    num_freq = round(Int64,max(50,2*f_range/(Δ_fP/4)))
    # num_freq = 1
    # num_freq = round(Int64, num_freq * pressure/100)
    # println("num_freq = ", num_freq)

    layer_unknown = n_rot*num_freq + n_vib
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
    # p_dist satisfies sum(p_dist) * df ~ 1.0 ;
    #p_dist = lorentz_dist(f_dist_ctr, Δ_f_NT, f_pump)
    fp_lasing = f_NT_ampl(f_dist_dir_lasing, Δ_fP, f_dir_lasing)
    fp_lasing = fp_lasing/sum(fp_lasing)
    fp_lasing = f_NT_normalized(f_dist_dir_lasing, Δ_fP, f_dir_lasing, df*f_dir_lasing/f₀)

    fp_ref_lasing = f_NT_ampl(f_dist_ref_lasing, Δ_fP, f_ref_lasing)
    fp_ref_lasing = fp_ref_lasing/sum(fp_ref_lasing)
    fp_ref_lasing = f_NT_normalized(f_dist_ref_lasing, Δ_fP, f_ref_lasing, df*f_ref_lasing/f₀)

    f_dirgain_dist = linspace(f_dir_lasing-40e6, f_dir_lasing+40e6, 500)
    f_refgain_dist = linspace(f_ref_lasing-40e6, f_ref_lasing+40e6, 500)
    dirgain = zeros(size(f_dirgain_dist))
    refgain = zeros(size(f_refgain_dist))
    # alpha_0 from eq (2.B.3) first line in unit m^-1:

    beta13 = 1.20 * sqrt(power)/radius

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

    totalNL0 = C4L * ntotal * f_G_0/2
    totalNU0 = C5U * ntotal * f_3_0/2
    # alpha_0 in m^-1
    alpha_0 = exp(-log(2)*((f_pump-f₀)/Δ_f₀D)^2)*sqrt(log(2)/pi)/Δ_f₀D *
              8*pi^3/3/h/c * (totalNL0 - totalNU0) *
              1e-36 * (0.2756^2*16.0/45) * f_pump * 1e-13
    alpha_r = alpha_0 * ones(num_layers)
    pump_IR = zeros(num_freq, num_layers)
    # in unit m^-3 microsec^-1
    pump_0 = 0
    # (1-exp(-alpha_0*L/100)) * power/(pi*(radius/100)^2*L/100 * h*f_pump *
    #          ntotal * f_G_0/2 * C4L) / norm_time
    pump_0 = 9.4e13 * power/(radius^2)/Δ_f₀D * (0.2756^2*16.0/45) *
                exp(-log(2)*((f_pump-f₀)/Δ_f₀D)^2)/norm_time

    Δr = radius/100 / num_layers # in m
    r_ext = linspace(0,radius/100,num_layers+1)
    r_int = 0.5*(r_ext[1:end-1] + r_ext[2:end]) # in m

    # kwall = WallRate(radius, pressure, r_int, ntotal, M, T, NA, v_avg, σ_GKC) + 1e-10
    kwall = zeros(num_layers)


    MFP = 0.732*T/pressure/σ_GKC # in cm
    # avg absolute vel in m/microsec (p5 in Henry's thesis)
    vel = v_avg/sqrt(2)/norm_time
    # diffusion coefficient in m^2/microsec. Einstein-Smoluchowski equation;
    D = 1/3 * vel * MFP * 1e-2 * D_factor

    kDD = 19.8 * pressure * σ_DD/ sqrt(T*M)

    J = [3,4,5,6,7,8,9,10,11]
    # in 1/microsec. rate_ij = rate_DD * prob. of ij collision:
    k98_G = kDD*Q_selectn_hl(J[9])/(1e3)
    k87_G = kDD*Q_selectn_hl(J[8])/(1e3)
    k76_G = kDD*Q_selectn_hl(J[7])/(1e3)
    k65_G = kDD*Q_selectn_hl(J[6])/(1e3)
    k54_G = kDD*Q_selectn_hl(J[5])/(1e3)
    k43_G = kDD*Q_selectn_hl(J[4])/(1e3)
    k32_G = kDD*Q_selectn_hl(J[3])/(1e3)
    k21_G = kDD*Q_selectn_hl(J[2])/(1e3)

    k89_G = kDD*Q_selectn_lh(J[8])/(1e3)
    k78_G = kDD*Q_selectn_lh(J[7])/(1e3)
    k67_G = kDD*Q_selectn_lh(J[6])/(1e3)
    k56_G = kDD*Q_selectn_lh(J[5])/(1e3)
    k45_G = kDD*Q_selectn_lh(J[4])/(1e3)
    k34_G = kDD*Q_selectn_lh(J[3])/(1e3)
    k23_G = kDD*Q_selectn_lh(J[2])/(1e3)
    k12_G = kDD*Q_selectn_lh(J[1])/(1e3)

    k98_3 = k98_G
    k87_3 = k87_G
    k76_3 = k76_G
    k65_3 = k65_G
    k54_3 = k54_G
    k43_3 = k43_G
    k32_3 = k32_G
    k21_3 = k21_G

    k89_3 = k89_G
    k78_3 = k78_G
    k67_3 = k67_G
    k56_3 = k56_G
    k45_3 = k45_G
    k34_3 = k34_G
    k23_3 = k23_G
    k12_3 = k12_G

    # K-swap rates -> goes to thermal pool, in 1/microsec
    k1a = 19.8*pressure*σ_SPT/sqrt(T*M)/(1e3)
    k2a = k1a
    k3a = k1a
    k4a = k1a
    k5a = k1a
    k6a = k1a
    k7a = k1a
    k8a = k1a
    k9a = k1a
    k10a = k1a
    k11a = k1a
    k12a = k1a
    k13a = k1a
    k14a = k1a
    k15a = k1a
    k16a = k1a
    k17a = k1a
    k18a = k1a

    return Params{DefaultT}(radius, pump_radius, L, L_eff, h, c, ev, kB, T, T_vA, T_vE, kBT, M, norm_time,
    σ_GKC, σ_DD, σ_SPT, σ_36, σ_VS,
    v_avg, kvs, #kvsplit,
    EG, E3, E6, E23, E36, E26, Q,
    f_G_0, f_3_0, f_6_0, f_GA, f_3A, f_6A, f_GE, f_3E, f_6E, f_23, f_36, f_26,
    C3L, C4L, C5L, C4U, C5U, g_L, g_U, NA,
    f₀, f_offset, f_pump, f_dir_lasing, f_ref_lasing, Δ_f₀D, f_range,
    n_rot, n_vib,
    mu0, eps0,
    mode_num, p_library, n0, t_spont, Δν_THz,
    pressure, power, powerF, powerB, averagePF, averagePB,
    num_layers, ntotal,
    k63A, k36A, k63E, k36E, netrate_36A, k3623, k2336, k2636, k3626, kro,
    Δ_fP, Δ_f_RabiF, Δ_f_RabiB, Δ_f_NTF, Δ_f_NTB,
    num_freq, layer_unknown, df, f_dist_end, f_dist_ctr, f_dist_ctrB,
    velocity, f_dist_dir_lasing, f_dist_ref_lasing,
    gauss_dist, SHBF, SHBB, fp_lasing, fp_ref_lasing,
    alpha_0, alpha_r, pump_0, pump_IR,
    Δr, r_int,
    kwall, MFP, D,
    k98_G, k87_G, k76_G, k65_G, k54_G, k43_G, k32_G, k21_G,
    k89_G, k78_G, k67_G, k56_G, k45_G, k34_G, k23_G, k12_G,
    k98_3, k87_3, k76_3, k65_3, k54_3, k43_3, k32_3, k21_3,
    k89_3, k78_3, k67_3, k56_3, k45_3, k34_3, k23_3, k12_3,
    k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a, k11a,
    k12a, k13a, k14a, k15a, k16a, k17a, k18a,
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

function Q_selectn_hl(J) #J->J-1
    K = 3
    return (J^2-K^2)/(J*(2*J+1))
end

function Q_selectn_lh(J) #J->J+1
    K = 3
    return ((J+1)^2-K^2)/((J+1)*(2*J+1))*exp(-ΔEr(J, J+1, 3, "V0", 35) * 1.6e-13)
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

function lorentz_dist(ν, Δ_f_NT, f_pump)
    return  1/π .* Δ_f_NT ./ ((ν - f_pump).^2 + Δ_f_NT^2)
end

function f_NT_ampl(ν, Δ_f_NT, f_pump)
    SHB = Δ_f_NT^2 ./ ((ν - f_pump).^2 + Δ_f_NT^2)
    # SHB = SHB/sum(SHB)
    return SHB
end

function f_NT_normalized(ν, Δ_f_NT, f_pump, df)
    SHB = 1/π * df * Δ_f_NT ./ ((ν - f_pump).^2 + Δ_f_NT^2)
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
  # spectrum = spectrum / max(spectrum)
  # spectrum = df * spectrum * p.beta13^2/γ^2
  return spectrum
end

function Qv(kB, T, script)
    if script == 1
      data = readdlm("/Users/fanwang/.julia/v0.4/OPFIR/src/E_vib.jl")
    else
      data = readdlm("/Users/fanwang/.julia/v0.4/OPFIR/src/E_vib.jl")
    end
    Q = 1.0
    for i in 1:size(data, 1)
        Q += data[i,2] * exp(-data[i, 1]/(kB*T*8065.73))
    end
    return Q
end
