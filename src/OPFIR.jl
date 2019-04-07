module OPFIR
using Distributions
using Cubature
using NLsolve, MUMPS

for f in ("parameters.jl", "parametersCO.jl", "CO.jl", "main.jl",
          "anderson_accel.jl", "fixedpoint.jl",
          "compute_rhs.jl", "compute_row_col_val.jl",
          "gain.jl", "pump.jl", "Tv_model.jl", "check.jl",
          "Sherman_Morrison.jl", "plotting.jl", "E_vib.jl")
    include(f)
end

end # module
