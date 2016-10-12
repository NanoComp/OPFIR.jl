using PyPlot

function plotting(p, sol)
  layers = [1, p.num_layers÷2, p.num_layers]
  for layer in layers
    # non-thermal population for rotational levels
    index_offset = p.layer_unknown*(layer-1)
    fig = figure()
    hold
    title("non-thermal population vs frequency")
    plot((p.f_dist_ref_lasing-p.f_ref_lasing)/1e6,
          sol[(1:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    plot((p.f_dist_ref_lasing-p.f_ref_lasing)/1e6,
          sol[(2:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    plot((p.f_dist_ref_lasing-p.f_ref_lasing)/1e6,
          sol[(3:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    plot((p.f_dist_dir_lasing-p.f_dir_lasing)/1e6,
          sol[(11:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    plot((p.f_dist_dir_lasing-p.f_dir_lasing)/1e6,
          sol[(12:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    plot((p.f_dist_dir_lasing-p.f_dir_lasing)/1e6,
          sol[(13:p.n_rot:p.n_rot*p.num_freq) + index_offset],"-")
    legend(["L-1", "L", "L+1", "U-1", "U", "U+1"])
    xlabel("frequency offset in lasing transition (MHz)")
    figname = string("nonthermal_pop_vs_freq_layer_",layer,".png")
    savefig(figname)

    fig = figure()
    hold
    ng0_J4 = (p.ntotal*p.f_G_0/2) * p.C4L .* p.gauss_dist
    plot(p.f_dist_ctr, ng0_J4)
    hold
    plot(p.f_dist_ctr, Nl_total_dist_layer(p, sol, layer))
    xlabel("frequency in pump transition")
    ylabel("population (m^-3)")
    legend(["no pump", "with pump"])
    title("Population of Ground state J=4, K=3")
    figname = string("popL_layer_",layer,".png")
    savefig(figname)

    fig = figure()
    hold
    n30_J5 = p.ntotal*p.f_3_0/2 * p.C5U .* p.gauss_dist
    n30_J4 = p.ntotal*p.f_3_0/2 * p.C4U .* p.gauss_dist
    plot(p.f_dist_ctr, n30_J5)
    hold
    plot(p.f_dist_ctr, Nl_total_dist_layer(p, sol, layer)) #
    xlabel("frequency in pump transition")
    ylabel("population (m^-3)")
    legend(["no pump", "with pump"])
    title("Population of V3 state J=5, K=3")
    figname = string("popU_layer_", layer, ".png")
    savefig(figname)

    fig = figure()
    hold
    plot((p.f_dist_dir_lasing-p.f_dir_lasing)/1e6, inv_U_dist_layer(p, sol, layer))
    xlabel("frequency offset in lasing transition (MHz)")
    title("pop inversion of U level")
    figname = string("invU_layer_", layer, ".png")
    savefig(figname)

    fig = figure()
    hold
    plot((p.f_dist_ref_lasing-p.f_ref_lasing)/1e6, inv_L_dist_layer(p, sol, layer))
    xlabel("frequency offset in lasing transition (MHz)")
    title("pop inversion of L level")
    figname = string("invL_layer_", layer, ".png")
    savefig(figname)

    fig = figure()
    hold
    plot((p.f_dist_ctr - p.f₀)/1e6, p.pump_IR[:, layer])
    xlabel("frequency offset in IR (MHz)")
    title("pump transition rate")
    figname = string("pumprate_", layer, ".png")
    savefig(figname)
  end

  fig = figure()
  hold
  plot(p.r_int, p.T_vA)
  plot(p.r_int, p.T_vE)
  xlabel("radial position")
  title("effective T")
  figname = string("TvA_TvE.png")
  savefig(figname)
end

function plot_gain_curve(centerfreq, freq, gain, filename)
  fig = figure()
  plot((freq-centerfreq)/1e6, gain)
  xlabel("offset frequency (MHz)")
  ylabel("gain (a.u.)")
  title(filename)
  savefig(string(filename, ".png"))
end
