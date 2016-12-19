
function update_alpha_from_N!(p, sol) # in m^-1
    tmp_factor = 8*π^3/3/p.h/p.c * 1e-36 * (0.2756^2*16.0/45) *
                  p.f_pump * 1e-13 / p.df
    for ri in 1:p.num_layers
        Nl_dist = Nl_total_dist_layer(p, sol, ri)
        Nu_dist = Nu_total_dist_layer(p, sol, ri)
        # abs_ability = f_NT_ampl(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump)
        abs_ability = f_NT_normalized(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump, p.df)
        # maxv, maxind = findmax(abs_ability)
        p.alpha_r[ri] = tmp_factor * sum((Nl_dist - p.g_L/p.g_U*Nu_dist).*abs_ability)
        # p.alpha_r[ri] = tmp_factor * (Nl_dist[maxind] - p.g_L/p.g_U*Nu_dist[maxind])
    end
    alpha = sum(p.alpha_r.*p.r_int)/sum(p.r_int)
    p.alpha_r = ones(p.num_layers) * alpha
    # println("alpha=", alpha)
end

function CavityAbsorption_FB(alpha_0, L, p)
    t1 = 0.95
    t2 = 0.96 * 0.95
    alphaL = alpha_0*L/100
    coeffF = (1 - exp(-alphaL)) / (1-t1*t2*exp(-2*alphaL))
    if p.MultRef == 0
      coeffF = 1-exp(-alphaL)
    end
    #coeffF = (1 - exp(-alphaL))
    coeffB = coeffF * exp(-alphaL)*t1
    return coeffF, coeffB
end

function AveragePower_FB(alpha_0, L, power)
    #coeffF, coeffB = CavityAbsorption_FB(alpha_0, L)
    power_local = power/(1-exp(-alpha_0*L/100))

    powerF = (1-exp(-alpha_0*L/100))/(alpha_0*L/100) * power_local
    # powerB = (1-exp(-alpha_0*L/100))/(alpha_0*L/100) * power_local
    powerB = powerF * exp(-alpha_0 * L/100)
    return powerF, powerB
end

function AveragePower2_FB(alpha_0, L, power)
    powerF = (1-exp(-alpha_0*L/100))/(alpha_0*L/100) * power
    powerB = powerF * exp(-alpha_0*L/100)
    return powerF, powerB
end

function update_Param_from_alpha!(p, sol)
    AbsorptionF = Array(Float64, p.num_layers)
    AbsorptionB = Array(Float64, p.num_layers)
    averagePF = Array(Float64, p.num_layers)
    averagePB = Array(Float64, p.num_layers)

    alpha_avg = sum(p.alpha_r .* p.r_int)/sum(p.r_int)
    L_eff = min(p.L, 200/alpha_avg) # in cm
    #L_eff = 15.0
    p.L_eff = L_eff

    for ri in 1:p.num_layers
      # percentage being absorbed in forward/backward direction
      # used for pump rate calculation
        AbsorptionF[ri], AbsorptionB[ri] = CavityAbsorption_FB(p.alpha_r[ri], p.L, p)
      # powerF/B for computing AC stark effect
        p.powerF[ri], p.powerB[ri] = AveragePower_FB(p.alpha_r[ri], p.L, p.power)

        p.averagePF[ri], p.averagePB[ri] = AveragePower2_FB(p.alpha_r[ri], p.L, p.power)

        p.Δ_f_RabiF[ri] = 0.38*sqrt(p.powerF[ri])/p.radius*1e6
        p.Δ_f_RabiB[ri] = 0.38*sqrt(p.powerB[ri])/p.radius*1e6
        p.Δ_f_RabiF[ri] = 0.38*sqrt(p.averagePF[ri])/p.radius*1e6
        p.Δ_f_RabiB[ri] = 0.38*sqrt(p.averagePB[ri])/p.radius*1e6

        # p.Δ_f_NTF[ri] = p.Δ_fP + p.Δ_f_RabiF[ri]
        # p.Δ_f_NTB[ri] = p.Δ_fP + p.Δ_f_RabiB[ri]
        p.Δ_f_NTF[ri] = sqrt(p.Δ_fP^2 + p.Δ_f_RabiF[ri]^2)
        p.Δ_f_NTB[ri] = sqrt(p.Δ_fP^2 + p.Δ_f_RabiB[ri]^2)
        if p.ACStark == 0
          p.Δ_f_NTF[ri] = p.Δ_fP
          p.Δ_f_NTB[ri] = p.Δ_fP
        end
        # p.SHBF[:, ri] = f_NT_ampl(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump)
        # fV = 0.5346*p.Δ_f_NTF[ri] + sqrt(0.2166*p.Δ_f_NTF[ri]^2 + p.Δ_f₀D^2)
        fV = p.Δ_f_NTF[ri]
        p.SHBF[:, ri] = f_NT_normalized(p.f_dist_ctr, fV, p.f_pump, p.df)
        # p.SHBB[:, ri] = f_NT_ampl(p.f_dist_ctrB, p.Δ_f_NTB[ri], p.f_pump)
        #fV = 0.5346*p.Δ_f_NTB[ri] + sqrt(0.2166*p.Δ_f_NTB[ri]^2 + p.Δ_f₀D^2)
        fV = p.Δ_f_NTB[ri]
        p.SHBB[:, ri] = f_NT_normalized(p.f_dist_ctrB, fV, p.f_pump, p.df)
    end

    for ri in 1:p.num_layers
        Nl_dist = Nl_total_dist_layer(p, sol, ri)
        Nu_dist = Nu_total_dist_layer(p, sol, ri)
        # abs_abilityF = f_NT_ampl(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump)
        # abs_abilityB = f_NT_ampl(p.f_dist_ctrB, p.Δ_f_NTB[ri], p.f_pump)
        abs_abilityF = f_NT_normalized(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump, p.df)
        abs_abilityB = f_NT_normalized(p.f_dist_ctrB, p.Δ_f_NTB[ri], p.f_pump, p.df)
        totalNLF = sum((Nl_dist-0*p.g_L/p.g_U*Nu_dist).*abs_abilityF)
        totalNLB = sum((Nl_dist-0*p.g_L/p.g_U*Nu_dist).*abs_abilityB)

        pumpregion = 1

        for vi in 1:p.num_freq
            if p.backward == 0
              p.pump_IR[vi, ri] =
              (p.SHBF[vi, ri] * AbsorptionF[ri] / totalNLF) *
              (1/pumpregion)^2 * p.power / (pi*(p.radius/100)^2*p.L/100 * p.h*p.f_pump) / p.norm_time
            elseif p.backward == 1
              p.pump_IR[vi, ri] =
              (p.SHBF[vi, ri] * AbsorptionF[ri] / totalNLF +
              p.SHBB[vi, ri] * AbsorptionB[ri] / totalNLB) *
              p.power / (pi*(p.radius/100)^2*p.L/100 * p.h*p.f_pump) / p.norm_time
            else
              error("p.backward can only be 0 or 1!")
            end

            # p.pump_IR[vi, ri] =
            # 8*π^3/(3*p.h^2*p.c) * 1e-36 * (0.2756^2*16.0/45) *p.power/(π*p.radius^2) * 1e-9 *
            # (p.SHBF[vi, ri]/p.df * p.powerF[ri]/p.power + p.SHBB[vi, ri]/p.df * p.powerB[ri]/p.power)/
            # p.norm_time
            p.pump_IR[vi, ri] =
            8*π^3/(3*p.h^2*p.c) * 1e-36 * (0.2756^2*16.0/45) *p.power/(π*p.radius^2) * 1e-9 *
            (p.SHBF[vi, ri]/p.df * p.averagePF[ri]/p.power +
             p.SHBB[vi, ri]/p.df * p.averagePB[ri]/p.power)/p.norm_time
            #            (p.SHBF[vi, ri]/p.df) / p.norm_time
            # assume constant pump rate all the time
            #p.pump_IR[vi, ri] = p.pump_0 * f_NT_ampl(p.f_dist_ctr[vi], p.Δ_f_NTF[ri], p.f_pump)

            if ri > p.num_layers * pumpregion
                p.pump_IR[vi, ri] = 0
            end
        end
    end

end
