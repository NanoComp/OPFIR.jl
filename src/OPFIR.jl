module OPFIR
using Distributions
using Cubature

for f in ("parameters.jl", "main.jl",
          "anderson_accel.jl", "fixedpoint.jl",
          "compute_rhs.jl", "compute_row_col_val.jl",
          "gain.jl")
    include(f)
end

end # module
