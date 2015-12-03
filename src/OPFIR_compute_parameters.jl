function OPFIR_compute_parameters(pressure, power, num_layers)
    ntotal::Real = 9.66e24 * pressure * 1e-3 / T; # with unit m^-3
    k63::Real = ntotal*v*σ_36*(1e-10)^2/norm_time/2; # in 1/microsec
    k36::Real = exp(-(E6-E3)/kBT) * k63 * 2; # in 1/microsec
    k3623::Real = ntotal*v*σ_36*(1e-10)^2/norm_time; # in 1/microsec
    k2336::Real = exp(-(E36-E23)/kBT) * k3623;
    k2636::Real = ntotal*v*σ_36*(1e-10)^2/norm_time;
    k3626::Real = exp(-(E26-E36)/kBT) * k2636;
    kro::Real = ntotal * v * σ_VS * (1e-10)^2 / norm_time; # in 1/microsec ;

    Δ_fP::Real = 15e6*(pressure/1e3); # Δ_fP = 15 MHz/Torr
    num_freq::Int64 = round(Int64,max(50,2*f_range/(Δ_fP/4)));
    # num_freq = 5;
    df::Real = 2.0 * f_range / num_freq;
    f_dist_end::Array = linspace(-f_range, f_range, num_freq + 1) + f₀;
    f_dist_ctr::Array = f_dist_end[1:end-1] + df/2;

    velocity::Array = (f_dist_ctr - f₀)/f₀; # in unit c
    f_dist_dir_lasing::Array = velocity * f_dir_lasing + f_dir_lasing;

    norm_dist = Normal(f₀, Δ_f₀D * sqrt(2*log(2)));
    pdf1 = pdf(norm_dist, f_dist_ctr);
    gauss_dist::Array = pdf1 / sum(pdf1);

    Δ_f_Rabi::Real = 0.45*sqrt(power)/radius*1e6; # HMHW in Hz
    #Δ_f_NT = sqrt(Δ_fP^2 + 4*Δ_f_Rabi^2);  # HMHW in Hz
    Δ_f_NT::Real = Δ_fP;

    p_dist = lorentz_dist(f_dist_ctr, Δ_f_NT, f_pump); # sum(p_dist) * df ~ 1.0 ;
    SHB::Array = f_NT_ampl(f_dist_ctr, Δ_f_NT, f_pump);
    # SHB = SHB/sum(SHB)
    # println("sum of SHB is ", sum(SHB))

    fp_lasing::Array = f_NT_ampl(f_dist_dir_lasing, Δ_f_NT, f_dir_lasing);
#     SHB = p_dist;

    # constants in pumping
    # in m-3 microsec-1 ;
    # pump0 = (3.487e-7/1e6/df) * power/(pi*(radius/100)^2)/(h*f₀)/norm_time / g_L;
    # pump0 = 9.4e13 * power/(radius^2)/Δ_f₀D * 0.098^2 * exp(-log(2)*((f_pump-f₀)/Δ_f₀D)^2)/norm_time
    pump0 = 9.4e13 * power/(radius^2)/Δ_f₀D * (0.2756^2*16.0/45) * exp(-log(2)*((f_pump-f₀)/Δ_f₀D)^2)/norm_time
    # println(pump0)
    # println((3.487e-7/1e6/df) * power/(pi*(radius/100)^2)/(h*f₀)/norm_time / g_L)
    pumpR::Array = pump0 * SHB;

    Δr::Real = radius/100 / num_layers; # in m
    r_ext = linspace(0,radius/100,num_layers+1);
    r_int::Array = 0.5*(r_ext[1:end-1] + r_ext[2:end]); # in m

    kwall::Array = WallRate(radius, pressure, r_int, ntotal) + 1e-10;

    MFP::Real = 0.732*T/pressure/σ_GKC; # in cm
    # avg absolute vel in m/microsec, p5 Henry's dissertation:
    vel = 1/3 * v/sqrt(2)/norm_time;
    # diffusion coefficient in m^2/microsec. Einstein-Smoluchowski equation;
    D::Real = vel * MFP * 1e-2;

    kDD::Real = 19.8 * pressure * σ_DD/ sqrt(T*M);

    J = [3,4,5,6,7,8,9,10,11];
    # in 1/microsec. rate_ij = rate_DD * prob. of ij collision:
    k98_G = kDD*Q_selectn_hl(J[9])/(1e3);
    k87_G = kDD*Q_selectn_hl(J[8])/(1e3);
    k76_G = kDD*Q_selectn_hl(J[7])/(1e3);
    k65_G = kDD*Q_selectn_hl(J[6])/(1e3);
    k54_G = kDD*Q_selectn_hl(J[5])/(1e3);
    k43_G = kDD*Q_selectn_hl(J[4])/(1e3);
    k32_G = kDD*Q_selectn_hl(J[3])/(1e3);
    k21_G = kDD*Q_selectn_hl(J[2])/(1e3);

    k89_G = kDD*Q_selectn_lh(J[8])/(1e3);
    k78_G = kDD*Q_selectn_lh(J[7])/(1e3);
    k67_G = kDD*Q_selectn_lh(J[6])/(1e3);
    k56_G = kDD*Q_selectn_lh(J[5])/(1e3);
    k45_G = kDD*Q_selectn_lh(J[4])/(1e3);
    k34_G = kDD*Q_selectn_lh(J[3])/(1e3);
    k23_G = kDD*Q_selectn_lh(J[2])/(1e3);
    k12_G = kDD*Q_selectn_lh(J[1])/(1e3);

    k98_3 = k98_G;
    k87_3 = k87_G;
    k76_3 = k76_G;
    k65_3 = k65_G;
    k54_3 = k54_G;
    k43_3 = k43_G;
    k32_3 = k32_G;
    k21_3 = k21_G;

    k89_3 = k89_G;
    k78_3 = k78_G;
    k67_3 = k67_G;
    k56_3 = k56_G;
    k45_3 = k45_G;
    k34_3 = k34_G;
    k23_3 = k23_G;
    k12_3 = k12_G;

    # Symmetry Preserving Thermalization (SPT) rates
    # K-swap rates -> goes to thermal pool
    k1a::Real = 19.8*pressure*σ_SPT/sqrt(T*M)/(1e3); # in 1/microsec
    k2a = k1a; # in 1/microsec
    k3a = k1a; # in 1/microsec
    k4a = k1a; # in 1/microsec
    k5a = k1a; # in 1/microsec
    k6a = k1a; # in 1/microsec
    k7a = k1a; # in 1/microsec
    k8a = k1a; # in 1/microsec
    k9a = k1a; # in 1/microsec
    k10a = k1a; # in 1/microsec
    k11a = k1a; # in 1/microsec
    k12a = k1a; # in 1/microsec
    k13a = k1a; # in 1/microsec
    k14a = k1a; # in 1/microsec
    k15a = k1a; # in 1/microsec
    k16a = k1a; # in 1/microsec
    k17a = k1a; # in 1/microsec
    k18a = k1a; # in 1/microsec ;

    return parameter(pressure, power, num_layers,
                     ntotal, k63, k36, k3623, k2336, k2636, k3626, kro,
                     Δ_fP, Δ_f_Rabi, Δ_f_NT,
                     num_freq, df,
                     f_dist_end, f_dist_ctr, velocity, f_dist_dir_lasing,
                     gauss_dist,
                     SHB, fp_lasing, pumpR,
                     Δr, r_int,
                     kwall,
                     MFP,
                     D,
                     k98_G, k87_G, k76_G, k65_G, k54_G, k43_G, k32_G, k21_G,
                     k89_G, k78_G, k67_G, k56_G, k45_G, k34_G, k23_G, k12_G,
                     k98_3, k87_3, k76_3, k65_3, k54_3, k43_3, k32_3, k21_3,
                     k89_3, k78_3, k67_3, k56_3, k45_3, k34_3, k23_3, k12_3,
                     k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a, k11a,
                     k12a, k13a, k14a, k15a, k16a, k17a, k18a,
                     1);
end
