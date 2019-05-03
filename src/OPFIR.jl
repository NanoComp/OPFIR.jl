module OPFIR

using Distributions
using Cubature
using NLsolve
using JLD
using ODE
using MUMPS

for f in ("parameters.jl", "parametersN2O.jl", "N2O.jl", "main.jl",
          "anderson_accel.jl", "fixedpoint.jl",
          "compute_rhs.jl", "compute_row_col_val.jl",
          "gain.jl", "pump.jl", "Tv_model.jl", "E_vib.jl",
          "utils.jl")
    include(f)
end

end # module
