function μij2_THz(p)
    if p.M == 44
        return 0.17^2/2
    end
    if p.M == 28
        return 0.12^2/2
    end
end

function μij2_IR(p)
    if p.M == 44 # n2o
        return 0.060373/2
    end
    if p.M == 28
        return 0.011568/2
    end
end

function Δfp(p)
    if p.M == 44
        return 4.0e6*(p.pressure/1e3)
    end
    if p.M == 28
        return 3.2e6*(p.pressure/1e3)
    end
end

function wallrate(u, Rcell, p)
    kw = p.MFP>p.radius ? 2u/3Rcell : 2/3*u*p.MFP/100/(p.radius/100)^2 # /(p.pressure/80)
    return kw
end

function kddrate(p)
    kdd = 19.8 * p.pressure * p.σ_DD/ sqrt(p.T*p.M)*1e3
end

function kSPTrate(p)
    kspt = p.ka[1]*1e6*0
end

function IRabsorption(αIR, Pqcl, p; model="TLM", multitrip=true)
    ## compute RHS (αIR) given αIR
    Lcell = p.L /100
    Rcell = p.radius/100
    tot_ir_abs = total_IR_abs(Pqcl, p, αIRgiven=true, αIR=αIR, multitrip=multitrip)
    Rpump = Pqcl * tot_ir_abs / (p.h*p.f₀) / (π*Rcell^2*Lcell)
    kw = wallrate(p.v_avg/sqrt(2), Rcell, p)
    kdd = kddrate(p)
    kspt = kSPTrate(p)
    if model == "TLM"
        N2 = Rpump/(2kdd+kw+kspt) *(kdd+kw+kspt)/(kw+kspt)
    elseif model == "FLM"
        N2 = Rpump/(1.5kdd+kw+kspt) *(kdd+kw+kspt)/(kw+kspt+0.5kdd)
    elseif model == "FLM2"
        N2 = Rpump * (kdd+kw+kspt)/((kdd+kw+kspt)^2 - 0.5kdd^2)
    end
    N1 = -Rpump/(kdd+kw)
    ΔN_IR0 = (p.CL * p.ntotal * p.f_G_0 - p.CU * p.ntotal * p.f_3_0)
    # S_nu0 = sqrt(log(2)/π)/p.Δ_f₀D
    S_nu = 1/π * 1/p.Δ_fP
    S_nu0 = 1/π * 1/sqrt(p.Δ_f₀D^2 + p.Δ_fP^2)
    α = 8*π^3/3/p.h/p.c * 1e-36 * μij2_IR(p) * p.f₀ * 1e-13 *
    (ΔN_IR0*S_nu0 + (N1-N2)*S_nu0) #
    # α = 8*π^3/3/p.h/p.c * 1e-36 * μij2_IR(p) * p.f₀ * 1e-13 *
    # (ΔN_IR0*sqrt(log(2)/π)/p.Δ_f₀D + (N1-N2)*sqrt(log(2)/π)/p.Δ_f₀D) #/π)
    return α
end

function compute_α(Pqcl, p; model="FLM", multitrip=true)
    # compute αIR self consistently
    α = nlsolve((fvec, x) -> begin
                        fvec[1] = IRabsorption(x[1], Pqcl, p, model=model, multitrip=multitrip) - x[1]
            end, [0.0], iterations=100)
    return α.zero[1]
end

function total_IR_abs(Pqcl, p; αIRgiven=false, αIR=0, model="FLM", multitrip=true)
    Lcell = p.L/100
    R1 = p.backmirrorR_IR
    R2 = 0.96*R1
    if !αIRgiven
        αIR = compute_α(Pqcl, p, model=model, multitrip=multitrip)
    end
    beta = exp(-αIR*Lcell)
    tot_ir_abs = (1-beta)*(1+R1*beta)/(1-R1*R2*beta^2)
#     tot_ir_abs = 1-beta^2
    return tot_ir_abs
end

function compute_th(p; given_α_cell=false, α_given=0.3, th0=0.1, Tpinhole=0.016, model="FLM")
    if given_α_cell
        α_cell = α_given
    else
        α_cell = OPFIR.ohmicloss(p)  - log(1-Tpinhole)/2/(p.L/100)
    end
    solvePth = nlsolve((fvec, x) -> begin
        # kw = wallrate(p.v_avg/sqrt(2), p.radius/100, p)
        #
        # kdd = kddrate(p) #19.8 * p.pressure * p.σ_DD/ sqrt(p.T*p.M)*1e3
        # kSPT = kSPTrate(p) # p.ka[1]*1e6 #* 0
        # muij2 = μij2_THz(p) * (3.3e-30)^2
        #
        fvec[1] = α_cell/gaincoeff(x[1], p, model=model) - 1.0
        # fvec[1] = α_cell/(2/(3*p.h^2*8.85e-12) * muij2/(p.v_avg*(p.radius/100)^2*p.L/100) *
        #     total_IR_abs(x[1], p, model=model) / p.f₀ / (tt*kdd+kw+kSPT)) - x[1]
            end, [th0], iterations=100)
    pth = solvePth.zero[1]
end

function gaincoeff(Pqcl, p; model="FLM")
    kw = wallrate(p.v_avg/sqrt(2), p.radius/100, p)
    kdd = kddrate(p) #19.8 * p.pressure * p.σ_DD/ sqrt(p.T*p.M)*1e3
    kspt = kSPTrate(p) #p.ka[1]*1e6 #* 0
    muij2 = μij2_THz(p) * (3.3e-30)^2
    if model == "TLM"
        tt = 2.0kdd+kw+kspt
    elseif model == "FLM"
        tt = 1.5kdd+kw+kspt
    elseif model == "FLM2"
        tt = ((kdd+kw+kspt)^2 - 0.5kdd^2)/(0.5kdd+kw+kspt)
    end

    Δν = p.Δ_f₀D*p.f_dir_lasing/p.f₀ + Δfp(p)
    λTHz = p.c/p.f_dir_lasing
    γ_0 = 2/(3*p.h^2*8.85e-12) * muij2/(λTHz*Δν*(p.radius/100)^2*p.L/100) *
            total_IR_abs(Pqcl, p, model=model) * Pqcl / p.f₀ / tt
    # γ_0 = 2/(3*p.h^2*8.85e-12) * muij2/(p.v_avg*(p.radius/100)^2*p.L/100) *
    #         total_IR_abs(Pqcl, p, model=model) * Pqcl / p.f₀ / (tt*kdd+kw+kSPT)
end

function PTHz(Pqcl, p; given_α_cell=false, α_given=0.3, Tpinhole = 0.016, model="TLM")
    νTHz = p.f_dir_lasing
    γ0 = gaincoeff(Pqcl, p, model=model)
    if given_α_cell
        α_cell = α_given
    else
        α_cell = OPFIR.ohmicloss(p)  - log(1-Tpinhole)/2/(p.L/100)
    end
    kdd = kddrate(p)
    kw = wallrate(p.v_avg/sqrt(2), p.radius/100, p)
    kspt = kSPTrate(p)

    tot_ir_abs = total_IR_abs(Pqcl, p, model=model)
    Rpump = Pqcl * tot_ir_abs / (p.h*p.f₀) / (π*(p.radius/100)^2*(p.L/100))
    if model == "TLM"
        ΔN = Rpump/(2kdd + kspt + kw)
        ks = kdd + 0.5kspt + 0.5kw
    elseif model == "FLM"
        ΔN = Rpump/(1.5kdd + kspt + kw)
        ks = 0.75kdd + 0.5kspt + 0.5kw
    elseif model == "FLM2"
        ΔN = Rpump*(0.5kdd+kw+kspt)/((kdd+kw+kspt)^2 - 0.5kdd^2)
        ks = ((kdd+kw+kspt)^2 - 0.5kdd^2)/(kdd+2(kw+kspt)-0.25kdd^2/(kdd+kw+kspt))
    end
    Δν = p.Δ_f₀D*p.f_dir_lasing/p.f₀ + Δfp(p)
    σ = 2π*μij2_THz(p)* (3.3e-30)^2/(3p.h*8.85e-12*p.c/νTHz*Δν)
    output_power= 0.5*(p.h*νTHz) * Tpinhole * (π*(p.radius/100)^2) *(γ0/α_cell-1) * ks/σ
end
