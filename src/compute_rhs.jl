function compute_rhs(rhs, p, sol)
    # pump_IR = Array(Float64, p.num_layers, p.num_freq)
    # update_alpha_from_N!(sol, p)
    # update_Param_from_alpha!(p)

    for ri in 1:p.num_layers
        for vi in 1:p.num_freq

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + 2
            rhs[row] = - p.pump_IR[vi, ri] * p.gauss_dist[vi] * p.ntotal *
                       (p.C4L * p.f_G_0/2 - p.C5U * p.f_3_0/2 * p.g_L/p.g_U)

            row = (ri-1)*p.layer_unknown + (vi-1)*p.n_rot + p.n_rot÷2 + 3
            rhs[row] = p.pump_IR[vi, ri] * p.gauss_dist[vi] * p.ntotal *
                       (p.C4L * p.f_G_0/2 - p.C5U * p.f_3_0/2 * p.g_L/p.g_U)
        end
    end

    if p.model_flag==2 && p.solstart_flag==1
        for ri in 1:p.num_layers
            # row of V3A level:
            row = (ri-1)*p.layer_unknown + p.num_freq*p.n_rot + 2
            # tot_pump = pump_total(p, sol, ri)
            # V3A -> VΣA net rate
            rhs[row] -= p.netrate_36A[ri]
            rhs[row+1] += p.netrate_36A[ri]
            # rhs[row] += tot_pump * p.f_6_0/(p.f_3_0+p.f_6_0)
            # rhs[row+1] += - tot_pump * p.f_6_0/(p.f_3_0+p.f_6_0)
        end
    end

end


function update_Param_from_alpha!(p)
    for ri in 1:p.num_layers
        DecayFactorF, DecayFactorB = CavityAbsorption_FB(p.alpha_r[ri], p.L)
        p.powerF[ri] = p.power * DecayFactorF
        p.powerB[ri] = p.power * DecayFactorB
        #Δ_f_NT = sqrt(Δ_fP^2 + Δ_f_Rabi^2)
        #Δ_f_NT = Δ_fP + Δ_f_Rabi
        #Δ_f_NT = Δ_fP # not include the Rabi oscillation
        p.Δ_f_RabiF[ri] = 0.45*sqrt(p.powerF[ri])/p.radius*1e6
        p.Δ_f_RabiB[ri] = 0.45*sqrt(p.powerB[ri])/p.radius*1e6
        p.Δ_f_NTF[ri] = p.Δ_fP + p.Δ_f_RabiF[ri]
        p.Δ_f_NTB[ri] = p.Δ_fP + p.Δ_f_RabiB[ri]

        p.SHBF[:, ri] = f_NT_ampl(p.f_dist_ctr, p.Δ_f_NTF[ri], p.f_pump)
        p.SHBB[:, ri] = f_NT_ampl(p.f_dist_ctrB, p.Δ_f_NTB[ri], p.f_pump)

        # pump rate in m-3 microsec-1:
        #pump0 = 9.4e13 * power/(radius^2)/Δ_f₀D * (0.2756^2*16.0/45) *
                # exp(-log(2)*((f_pump-f₀)/Δ_f₀D)^2)/norm_time
        #pumpR = pump0 * SHB
        # pump_const_term = 9.4e13 * p.power/(p.radius^2)/p.Δ_f₀D * (0.2756^2*16.0/45) *
                    # exp(-log(2)*((p.f_pump-p.f₀)/p.Δ_f₀D)^2)/p.norm_time
        # p.pumpR = pump_const_term * (DecayFactorF * p.SHBF + DecayFactorB * p.SHBB)
    end

    for ri in 1:p.num_layers
        for vi in 1:p.num_freq
            p.pump_IR[vi, ri] = p.pump_0/p.power *
                    (p.powerF[ri] * p.SHBF[vi, ri] + p.powerB[ri] * p.SHBB[vi, ri])
        end
    end

end

function update_alpha_from_N!(sol, p)
    nonthermalNL = 0.0
    nonthermalNU = 0.0
    thermalNL = 0.0
    thermalNU = 0.0
    totalNL = nonthermalNL + thermalNL
    totalNL0 = p.C4L * p.ntotal * p.f_G_0/2
    totalNU = nonthermalNU + thermalNU
    totalNU0 = p.C5U * p.ntotal * p.f_3_0/2

    for ri in 1:p.num_layers
        nonthermalNL_dist = Nl_NT_dist_layer(p, sol, ri)
        nonthermalNL = sum(nonthermalNL_dist)

        thermalNL_dist = Nl_T_dist_layer(p, sol, ri)
        thermalNL = sum(thermalNL_dist)

        nonthermalNU_dist = Nu_NT_dist_layer(p, sol, ri)
        nonthermalNU = sum(nonthermalNU_dist)

        thermalNU_dist = Nu_T_dist_layer(p, sol, ri)
        thermalNU = sum(thermalNU_dist)

        totalNL = nonthermalNL + thermalNL
        totalNU = nonthermalNU + thermalNU
        p.alpha_r[ri] = (totalNL-totalNU) / (totalNL0-totalNU0) * p.alpha_0
    end

end


function pumpR_SHB!(sol, p)
    alpha_0 = p.alpha_0;

end

function pump_total(p, sol, layer)
    L_tot = Nl_total_dist_layer(p, sol, layer)
    U_tot = Nu_total_dist_layer(p, sol, layer)
    pump_dist = p.pumpR .* (L_tot-p.g_L/p.g_U*U_tot)
    return  sum(pump_dist)
end
