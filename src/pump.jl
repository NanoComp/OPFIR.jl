# function Nlu_total_layer(p, sol, ri)
#     nonthermalNL_dist = Nl_NT_dist_layer(p, sol, ri)
#     nonthermalNL = sum(nonthermalNL_dist)
#
#     thermalNL_dist = Nl_T_dist_layer(p, sol, ri)
#     thermalNL = sum(thermalNL_dist)
#
#     nonthermalNU_dist = Nu_NT_dist_layer(p, sol, ri)
#     nonthermalNU = sum(nonthermalNU_dist)
#
#     thermalNU_dist = Nu_T_dist_layer(p, sol, ri)
#     thermalNU = sum(thermalNU_dist)
#
#     totalNL = nonthermalNL + thermalNL
#     totalNU = nonthermalNU + thermalNU
#     return totalNL, totalNU
# end

function update_alpha_from_N!(p, sol)
    tmp_factor = 8*π^3/3/p.h/p.c * 1e-36 * (0.2756^2*16.0/45) * p.f_pump * 1e-13 / p.df
    # totalNL0 = p.C4L * p.ntotal * p.f_G_0/2
    # totalNU0 = p.C5U * p.ntotal * p.f_3_0/2
    for ri in 1:p.num_layers
        Nl_dist = Nl_total_dist_layer(p, sol, ri)
        Nu_dist = Nu_total_dist_layer(p, sol, ri)
        # abs_ability = f_NT_ampl(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump)
        abs_ability = f_NT_normalized(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump, p.df)
        # maxv, maxind = findmax(abs_ability)
        p.alpha_r[ri] = tmp_factor * sum((Nl_dist - p.g_L/p.g_U*Nu_dist).*abs_ability)
        # p.alpha_r[ri] = tmp_factor * (Nl_dist[maxind] - p.g_L/p.g_U*Nu_dist[maxind])
    end
    println(p.alpha_r)
end

function CavityAbsorption_FB(alpha_0, L)
    t1 = 0.95
    t2 = 0.96 * 0.95
    alphaL = alpha_0*L/100
    coeffF = (1 - exp(-alphaL)) / (1-t1*t2*exp(-2*alphaL))
    coeffB = coeffF * exp(-alphaL)*t1
    return coeffF, coeffB
end

function AveragePower_FB(alpha_0, L, power)
    coeffF, coeffB = CavityAbsorption_FB(alpha_0, L)
    coeffF = coeffF/(alpha_0*L/100) * power
    coeffB = coeffB/(alpha_0*L/100) * power
    return coeffF, coeffB
end

function update_Param_from_alpha!(p, sol)
    AbsorptionF = Array(Float64, p.num_layers)
    AbsorptionB = Array(Float64, p.num_layers)

    alpha_avg = sum(p.alpha_r .* p.r_int)/sum(p.r_int)
    L_eff = min(p.L, 200/alpha_avg)
    #L_eff = 15.0
    p.L_eff = L_eff
    println("L_eff = ", L_eff)
    
    for ri in 1:p.num_layers
        AbsorptionF[ri], AbsorptionB[ri] = CavityAbsorption_FB(p.alpha_r[ri], L_eff)
        p.powerF[ri], p.powerB[ri] = AveragePower_FB(p.alpha_r[ri], L_eff, p.power)
        # p.powerF[ri] = p.power
        # p.powerB[ri] = p.power
        p.Δ_f_RabiF[ri] = 0.45*sqrt(p.powerF[ri])/p.radius*1e6
        p.Δ_f_RabiB[ri] = 0.45*sqrt(p.powerB[ri])/p.radius*1e6
        p.Δ_f_NTF[ri] = p.Δ_fP + p.Δ_f_RabiF[ri]
        p.Δ_f_NTB[ri] = p.Δ_fP + p.Δ_f_RabiB[ri]

        # p.Δ_f_NTF[ri] = p.Δ_fP
        # p.Δ_f_NTB[ri] = p.Δ_fP

        # p.SHBF[:, ri] = f_NT_ampl(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump)
        p.SHBF[:, ri] = f_NT_normalized(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump, p.df)
        # p.SHBB[:, ri] = f_NT_ampl(p.f_dist_ctrB, p.Δ_f_NTB[ri], p.f_pump)
        p.SHBB[:, ri] = f_NT_normalized(p.f_dist_ctrB, p.Δ_f_NTB[ri], p.f_pump, p.df)
    end

    for ri in 1:p.num_layers
        Nl_dist = Nl_total_dist_layer(p, sol, ri)
        # abs_abilityF = f_NT_ampl(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump)
        # abs_abilityB = f_NT_ampl(p.f_dist_ctrB, p.Δ_f_NTB[ri], p.f_pump)
        abs_abilityF = f_NT_normalized(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump, p.df)
        abs_abilityB = f_NT_normalized(p.f_dist_ctrB, p.Δ_f_NTB[ri], p.f_pump, p.df)
        totalNLF = sum(Nl_dist.*abs_abilityF)
        totalNLB = sum(Nl_dist.*abs_abilityB)

        pumpregion = 1
        
        for vi in 1:p.num_freq
            p.pump_IR[vi, ri] =
            (p.SHBF[vi, ri] * AbsorptionF[ri] / totalNLF) *
            (1/pumpregion)^2 * p.power / (pi*(p.radius/100)^2*L_eff/100 * p.h*p.f_pump) / p.norm_time
            # p.pump_IR[vi, ri] =
            # (p.SHBF[vi, ri] * AbsorptionF[ri] / totalNLF +
            # p.SHBB[vi, ri] * AbsorptionB[ri] / totalNLB) *
            # p.power / (pi*(p.radius/100)^2*p.L/100 * p.h*p.f_pump) / p.norm_time
            # 9.4e13 * p.power/(p.radius^2)/p.Δ_f₀D * (0.2756^2*16.0/45) *
            #             exp(-log(2)*((p.f_pump-p.f₀)/p.Δ_f₀D)^2)/p.norm_time
            # p.pump_IR[vi, ri] *= p.SHBF[vi, ri]
            #

            if ri > p.num_layers * pumpregion
                p.pump_IR[vi, ri] = 0
            end
        end
    end

end
