# TODO:

type Params{T<:Real}
    radius::T # radius in cm
    L::T # length in cm
end

function Params(T=Float64; radius=0.25, L=15)
    return Params{T}(radius, L)
end
