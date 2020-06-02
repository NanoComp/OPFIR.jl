using LinearAlgebra

"""
    andersonaccel!(g!, x; m, norm=Base.norm, reltol=sqrt(ɛ), abstol=0, maxiters)

Solve the fixed-point problem \$g(x)=x\$ to given relative and absolute tolerances,
returning the solution \$x\$, via Anderson acceleration of a fixed-point
iteration starting at `x` (which must be of the correct type to hold the result).

This implementation operates in-place as much as possible.   You supply
a function `g!(y,x)` that computes \$y = g(x)\$, and `x` is overwritten
by the solution.

`m` is the number of iterates that are "remembered" by the Anderson acceleration.
`m=1` corresponds to unaccelerated fixed-point iteration.  The default is `m=10`
or `length(x) ÷ 2`, whichever is smaller.

The iteration halts when `norm(Δx) ≤ reltol*norm(x) + abstol` or `maxiters` iterations
is reached.  By default, `reltol` is the square root of the precision of `x`,
`abstol` is zero, and `maxiters` is `typemax(Int)`.  `norm` defaults to the built-in
`norm(x)` function (the L₂ norm), but can be changed via the `norm` keyword.
"""
function andersonaccel!(g!, x::Union{AbstractVector{R},AbstractVector{Complex{R}}};
                                          m = max(1, min(10, (length(x)+1)÷2)), # max #steps to remember
                                          norm=LinearAlgebra.norm, # norm to use for stop tolerance
                                          reltol::Real = sqrt(eps(real(R))),
                                          abstol::Real=0, maxiter::Int=typemax(Int)) where R <: AbstractFloat
    m < 1 && throw(ArgumentError("m=$m < 1 is not allowed"))
    reltol < 0 && throw(ArgumentError("reltol=$reltol < 0 is not allowed"))
    abstol < 0 && throw(ArgumentError("abstol=$abstol < 0 is not allowed"))
    n = length(x)
    m > max(1,n+1) && throw(ArgumentError("m=$m > n-1 = $n-1 is not allowed"))

    T = eltype(x)
    y = Array{T}(undef, n)

    if m == 1 # simple fixed-point iteration, no memory
        for k = 1:maxiter
            g!(y, x)
            for i = 1:n
                xnew = y[i]
                y[i] = xnew - x[i]
                x[i] = xnew
            end
            if norm(y) <= reltol*norm(x) + abstol
                return x
            end
        end
        return x
    end

    # pre-allocate all of the arrays we will need.  The
    # goal is to allocate once and re-use the storage
    # during the iteration by operating in-place.
    f = Array{T}(undef, n)
    X = Array{T}(undef, n, m-1)
    F = Array{T}(undef, n, m-1)
    Q = Array{T}(undef, n, m-1) # space for QR factorization
    γ = Array{T}(undef, max(n,m-1)) # not m-1, to store rhs (f) and overwrite in-place via A_ldiv_B!

    # first iteration is just x₂ = g(x₁) = y₁
    g!(y, x)
    kcol = 1
    for i = 1:n
        X[i,kcol] = f[i] = y[i] - x[i]
        x[i] = y[i]
    end
    if norm(f) <= reltol*norm(x) + abstol
        return x
    end

    for k = 2:maxiter
        g!(y, x)
        y[1]==0 && return 0
        for i = 1:n
            f_old = f[i]
            γ[i] = f[i] = y[i] - x[i]
            F[i,kcol] = f[i] - f_old
        end

        # construct subarrays to work in-place on a
        # subset of the columns
        mₖ = min(m, k)
        Xₖ = view(X, 1:n, 1:mₖ-1)
        Fₖ = view(F, 1:n, 1:mₖ-1)
        γₖ = view(γ, 1:mₖ-1)

        # use this once Julia issue #13728 is fixed:
        # Qₖ = view(Q, 1:n, 1:mₖ-1)
        # QR = qrfact!(copy!(Qₖ, Fₖ), Val{true})
        QR = m == mₖ ? qr!(copyto!(Q, Fₖ), Val(true)) : qr(Fₖ, Val(true))
        # LinearAlgebra.A_ldiv_B!(QR, γ) # overwrites γₖ in-place with Fₖ \ f
        LinearAlgebra.ldiv!(QR, γ)
        # We replace columns of F and X with the new
        # data in-place.  Rather than always appending
        # the new data in the last column, we cycle
        # through the m-1 columns periodically.
        kcol = (kcol % (m-1)) + 1 # next column of F, X to update

        # x = x + f - (Xₖ + Fₖ)*γₖ, updating in-place
        # (also update X[:,kcol+1] with new Δx, and set y = Δx)
        for i = 1:n
            xnew = y[i] # == x[i] + f[i]
            for j = 1:mₖ-1
                xnew -= (Xₖ[i,j] + Fₖ[i,j])*γₖ[j]
            end
            X[i,kcol] = y[i] = xnew - x[i]
            x[i] = xnew
        end
        if norm(y) <= reltol*norm(x) + abstol
            return x
        end
    end

    return x
end

"""
    andersonaccel(g, x; m, norm=LinearAlgebra.norm, reltol=sqrt(ɛ), abstol=0, maxiters)

Solve the fixed-point problem \$g(x)=x\$ to given relative and absolute tolerances,
returning the solution \$x\$, via Anderson acceleration of a fixed-point
iteration starting at `x` and given the function `g(x)`.

The keyword parameters are the same as for `andersonaccel`: they specify
the "memory" `m` of the algorithm, the relative (`reltol`) and absolute
(`abstol`) stopping tolerances in the given `norm`, and the maximum
number of itertions (`maxiters`).
"""
andersonaccel(g, x::AbstractVector{T}; kws...) where T <: Number =
    andersonaccel!((y,x) -> copyto!(y, g(x)),
                   copyto!(Array{typeof(float(one(T)))}(undef,length(x)), x);
                   kws...)

andersonaccel(g, x::T; kws...) where T <: Number = andersonaccel(x -> g(x[1]), [x]; kws...)[1]
