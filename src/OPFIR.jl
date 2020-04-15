module OPFIR

using Distributions
using Cubature
using NLsolve
using JLD
using ODE
using SparseArrays, LinearAlgebra, StatsBase, SpecialFunctions

for f in ("ParamsLinearMolecule.jl", "LinearMolecule.jl",
          "main.jl", "fixedpoint.jl",
          "anderson_accel.jl",
          "compute_rhs.jl", "compute_row_col_val.jl",
          "gain.jl", "pump.jl", "E_vib.jl",
          "three_level_model.jl")
    include(f)
end

end # module
