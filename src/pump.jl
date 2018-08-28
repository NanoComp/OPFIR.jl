function alphaback(p, sol) # in m^-1
    tmp_factor = 8*π^3/3/p.h/p.c * 1e-36 * dipolemom(p) * p.f_pump * 1e-13 / p.df

    alphab = zeros(p.num_layers)
    for ri in 1:p.num_layers
        Nl_dist = Nl_total_dist_layer(p, sol, ri)
        Nu_dist = Nu_total_dist_layer(p, sol, ri)
        abs_ability = f_NT_normalized(p.f_dist_ctrB, p.Δ_f_NTB[ri], p.f_pump, p.df)
        alphab[ri] = tmp_factor * sum((Nl_dist - p.g_L/p.g_U*Nu_dist).*abs_ability)
    end
    return sum(alphab.*p.r_int)/sum(p.r_int)*ones(p.num_layers)
end

function AveragePower2_FB(alphaf, alphab, L, power)
    powerF = (1-exp(-alphaf*L/100))/(alphaf*L/100) * power
    powerB = (1-exp(-alphab*L/100))/(alphab*L/100) * power * exp(-alphaf*L/100)
    return powerF, powerB
end

function comppumprate!(p, alphab)
    # numtrips = ceil(Int, log(0.01)/log(exp(-p.alpha_r[1]*p.L/100)*exp(-alphab[1]*p.L/100) * 0.95^2*0.96))
    R1 = 0.95
    R2 = 0.95 * 0.96
    triploss = exp(-p.alpha_r[1]*p.L/100)*exp(-alphab[1]*p.L/100) * R1*R2
    alphaforw = p.alpha_r[1]
    alphaback = alphab[1]
    pavgforw = (1-exp(-alphaforw*p.L/100))/(alphaforw*p.L/100) * p.power / (1-triploss)
    pavgback = (1-exp(-alphaback*p.L/100))/(alphaback*p.L/100) * p.power / (1-triploss) * exp(-alphaforw*p.L/100) * R1

    for ri in 1:p.num_layers
        p.Δ_f_RabiF[ri] = 0.38*sqrt(pavgforw)/p.radius*1e6
        p.Δ_f_RabiB[ri] = 0.38*sqrt(pavgback)/p.radius*1e6
        p.Δ_f_NTF[ri] = p.Δ_fP + p.Δ_f_RabiF[ri]
        p.Δ_f_NTB[ri] = p.Δ_fP + p.Δ_f_RabiB[ri]
        fV = p.Δ_f_NTF[ri]
        p.SHBF[:, ri] = f_NT_normalized(p.f_dist_ctr, fV, p.f_pump, p.df)
        fV = p.Δ_f_NTB[ri]
        p.SHBB[:, ri] = f_NT_normalized(p.f_dist_ctrB, fV, p.f_pump, p.df)

        p.pump_IR[:, ri] =
        8*π^3/(3*p.h^2*p.c)*1e-36*dipolemom(p)/(π*p.radius^2)*1e-9 *
        (p.SHBF[:, ri]/p.df * pavgforw +
         p.SHBB[:, ri]/p.df * pavgback)/p.norm_time
     end
end

function update_alpha_from_N!(p, sol) # in m^-1
    tmp_factor = 8*π^3/3/p.h/p.c * 1e-36 * dipolemom(p) * p.f_pump * 1e-13 / p.df

    for ri in 1:p.num_layers
        Nl_dist = Nl_total_dist_layer(p, sol, ri)
        Nu_dist = Nu_total_dist_layer(p, sol, ri)
        abs_ability = f_NT_normalized(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump, p.df)
        p.alpha_r[ri] = tmp_factor * sum((Nl_dist - p.g_L/p.g_U*Nu_dist).*abs_ability)
    end
    alpha = sum(p.alpha_r.*p.r_int)/sum(p.r_int)
    p.alpha_r = ones(p.num_layers) * alpha
    # println("alpha=", alpha)
end

function update_Param_from_alpha!(p, sol)
    alpha_avg = sum(p.alpha_r .* p.r_int)/sum(p.r_int)
    L_eff = min(p.L, 200/alpha_avg) # in cm
    #L_eff = 15.0
    p.L_eff = L_eff

    alphab = sum(alphaback(p, sol).*p.r_int)/sum(p.r_int) * ones(p.num_layers)
    comppumprate!(p, alphab)
end

function dipolemom(p)
    matele = 0.
    if p.pumpbranch == "R"
        matele = 0.060373*((p.JL+1)^2-p.K0^2)/(p.JL+1)/(2p.JL+1)
    elseif p.pumpbranch == "Q"
        matele = 0.060373*p.K0^2/(p.JL+1)/p.JL
    else
        matele = 0.060373*(p.JL^2-p.K0^2)/p.JL/(2p.JL+1)
    end
    return matele
end
