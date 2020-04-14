mutable struct ParamsLinearMolecule3{T<:Real}
    ## cavity property
    radius::T   # radius in cm
    pump_radius::T
    L::T        # length in cm
    L_eff::T    # effective length in cm
    T::T        # cavity temperature in K
    cavitywall::AbstractString # material of the cavity wall, such as "Cu", "Au", etc
    cavitymode::AbstractString
    frontmirrorT_THz::T # front mirror transmission coeff
    frontmirrorT_IR::T  # front mirror transmission coeff for IR
    backmirrorR_IR::T  # back mirror reflection for IR
    cavityloss::T # in m-1
    lossfactor::T
    lasinglevel::Union{AbstractString, Char}

    ## physical constants
    h::T        # Planck constant
    c::T        # speed of light in m/s
    ev::T       # in J
    kB::T       # in eV/K
    norm_time::T# normalization to micro second
    NA::T
    mu0::T
    eps0::T

    kBT::T      # in cm^-1

    ## molecule properties
    name::AbstractString # molecule name, such as "N2O", "CO", "HCN"
    M::T        # molecular mass in AMU
    σ_GKC::T    # in Angstrom^2, gas kinetic collision cross section
    σ_DD::T     # dipole - dipole
    σ_SPT::T

    pressure::T
    ntotal::T

    kDD::T
    kDDmat::AbstractArray{T,2}
    ka::T

    pbroadenwidth::T # per Torr, for example, 4.0e6 Hz for N2O
    pbfactor::T
    Δ_fP0::T
    ΔfQCL::T
    Δ_fP::T

    beta13::T
    Δ_f_RabiF::AbstractVector{T}
    Δ_f_RabiB::AbstractVector{T}
    Δ_f_NTF::AbstractVector{T}
    Δ_f_NTB::AbstractVector{T}


    v_avg::T    # in m/sec, average relative velocity between molecules
    MFP::T
    D::T
    D_factor::T

    B::T
    DJ::T
    JL::Integer #pump: JL in V0 ->JU in V3
    JU::Integer
    K0::Integer # = 0
    pumpbranch::AbstractString #'R' or 'P' or 'Q'
    CL::T
    CL1::T
    CU::T
    CU1::T
    g_L::T
    g_U::T

    f₀::T
    Δ_f₀D::T
    f_offset::T
    f_pump::T

    dipolemoment::T     # for spontaneous emission
    dipolederivative::T # for IR absorption


    alpha_0::T
    alpha_r::AbstractVector{T}
    # vibrational levels
    EG::T
    E3::T
    Q::T        # partition function
    f_G_0::T
    f_3_0::T

    f_dir_lasing::T
    f_ref_lasing::T

    t_spont::T

    Js::AbstractVector{T}
    n_rot::Integer
    n_vib::Integer
    # rotpopfr::AbstractVector{T}
    cj::AbstractVector{T}
    frel_vib::AbstractVector{T} # fraction of vibrational levels
    kVVmat::AbstractArray{T,2}  # VV collision constant matrix

    n0::T #refractive index

    ### power related:
    power::T        # in unit W
    powerF::AbstractVector{T}
    powerB::AbstractVector{T}
    averagePF::AbstractVector{T}
    averagePB::AbstractVector{T}
    pump_0::T
    pump_IR::AbstractArray{T,2}

    ## numerical (discretization) parameters
    f_range::T
    velocity::AbstractVector{T}
    f_dist_end::AbstractVector{T}
    f_dist_ctr::AbstractVector{T}
    f_dist_ctrB::AbstractVector{T}
    df::T
    num_freq::Integer
    layer_unknown::Integer

    f_dist_dir_lasing::AbstractVector{T}
    f_dist_ref_lasing::AbstractVector{T}
    gauss_dist::AbstractVector{T}
    SHBF::AbstractArray{T,2}
    SHBB::AbstractArray{T,2}

    Δr::T
    r_int::AbstractVector{T}
    num_layers::Integer

    solstart_flag::Integer
    mumps_solver::Integer

    # Wi::T
    WiU::T
    WiL::T

    ### important output quantities
    taus::T
    THzpower::T
    gaincoeff::T
end
