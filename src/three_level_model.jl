function μij2()
    return 0.060373/2
end

function wallrate(u, Rcell, p)
    # return kw = p.v_avg * p.MFP/100/(p.radius/100)^2
    return 2u/3Rcell
end

function kddrate(p)
    kdd = 19.8 * p.pressure * p.σ_DD/ sqrt(p.T*p.M)*1e3
end

function kSPTrate(p)
    kspt = p.ka[1]*1e6
end

function IRabsorption(αIR, Pqcl, p)
    Lcell = p.L /100
    Rcell = p.radius/100
    tot_ir_abs = total_IR_abs(Pqcl, p, αIRgiven=true, αIR=αIR)
    Rpump = Pqcl * tot_ir_abs / (p.h*p.f₀) / (π*Rcell^2*Lcell)
    kw = wallrate(p.v_avg, Rcell, p)
    N2pN3 = Rpump/kw
    ΔN_IR0 = (p.CL * p.ntotal * p.f_G_0 - p.CU * p.ntotal * p.f_3_0)
    α = 8*π^3/3/p.h/p.c * 1e-36 * μij2() * p.f₀ * 1e-13 *
    (ΔN_IR0*sqrt(log(2)/π)/p.Δ_f₀D - N2pN3*sqrt(log(2)/π)/p.Δ_f₀D) #/π)
    return α
end

function compute_α(Pqcl, p)
    α = nlsolve((fvec, x) -> begin
                        fvec[1] = IRabsorption(x[1], Pqcl, p) - x[1]
            end, [0.0], iterations=100)
    return α.zero[1]
end

function total_IR_abs(Pqcl, p; αIRgiven=false, αIR=0)
    Lcell = p.L/100
    R1 = 0.95
    R2 = 0.96*R1
    if !αIRgiven
        αIR = compute_α(Pqcl, p)
    end
    beta = expm(-αIR*Lcell)
    tot_ir_abs = (1-beta)*(1+R1*beta)/(1-R1*R2*beta^2)
#     tot_ir_abs = 1-beta^2
    return tot_ir_abs
end

function compute_th(p; α_cell=0.3)
    solvePth = nlsolve((fvec, x) -> begin
        kw = wallrate(p.v_avg, p.radius/100, p)

        kdd = kddrate(p) #19.8 * p.pressure * p.σ_DD/ sqrt(p.T*p.M)*1e3
        kSPT = kSPTrate(p) # p.ka[1]*1e6 #* 0
        muij2 = μij2() * (3.3e-30)^2
        fvec[1] = α_cell/(2/(3*p.h^2*8.85e-12) * muij2/(p.v_avg*(p.radius/100)^2*p.L/100) *
            total_IR_abs(x[1], p) / p.f₀ / (2kdd+kw+kSPT)) - x[1]
            end, [0.], iterations=100)
    pth = solvePth.zero[1]
end

function gaincoeff(Pqcl, p)
    kw = wallrate(p.v_avg, p.radius/100, p)
    kdd = kddrate(p) #19.8 * p.pressure * p.σ_DD/ sqrt(p.T*p.M)*1e3
    kSPT = kSPTrate(p) #p.ka[1]*1e6 #* 0
    muij2 = μij2() * (3.3e-30)^2
    γ_0 = 2/(3*p.h^2*8.85e-12) * muij2/(p.v_avg*(p.radius/100)^2*p.L/100) *
            total_IR_abs(Pqcl, p) * Pqcl / p.f₀ / (2kdd+kw+kSPT)
end

function PTHz(Pqcl, p; α_cell=0.3)
    νTHz = p.f_dir_lasing
    Tpinhole = 0.016
    γ0 = gaincoeff(Pqcl, p)
#     println(γ0)
    kdd = kddrate(p)
    kw = wallrate(p.v_avg, p.radius/100, p)
    kspt = kSPTrate(p)
    ks = kdd + kspt + 0.5kw
    tot_ir_abs = total_IR_abs(Pqcl, p)
    Rpump = Pqcl * tot_ir_abs / (p.h*p.f₀) / (π*(p.radius/100)^2*(p.L/100))
    ΔN = Rpump/(2kdd + 2kspt + kw)
    output_power= 0.5*(p.h*νTHz) * Tpinhole * (π*(p.radius/100)^2) *(γ0/α_cell-1) * ks*ΔN/γ0
end
