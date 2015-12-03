module OPFIR

include("parameters.jl")
for f in ("Lorentz_distribution.jl", "Sherman_Morrison.jl",
          "OPFIR_compute_parameters.jl", "anderson_accel.jl",
          "OPFIR_compute_rhs.jl", "population_inv.jl",
          "OPFIR_compute_row_col_val.jl", "post_gain.jl",
          "OPFIR_fixedpoint.jl", "post_plotting.jl",
          "OPFIR_func.jl", "wallrates.jl", "Q_select.jl")
    include(f)
end

end
