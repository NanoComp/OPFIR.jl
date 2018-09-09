type ParamsN2O{T<:Real}
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

    rotpopfr::AbstractVector
    cj::AbstractVector
end
