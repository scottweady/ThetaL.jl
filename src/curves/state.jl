# The θ–L state of a closed curve (`ThetaLState`: φ = θ − α, L, x₀) and what follows from it: the curve as a
# trigonometric interpolant, its tangent angles and curvature, the normal and tangential parts of a velocity, the
# area, centroid and length, and how well N resolves it. (Loops on concrete vectors: compile time, src/grid/fourier.jl.)

# What follows from a state and is used more than once (the half spectrum of φ, θ_α = 1 + φ_α, the grid points),
# computed at most once per state on first use, or filled in by the step that made the state. Not part of its
# value: a state is never changed after it is made.
mutable struct _StateCache
    φh::Union{Nothing,Vector{ComplexF64}}
    θα::Union{Nothing,Vector{Float64}}
    points::Union{Nothing,Vector{Vec2}}
end
_StateCache() = _StateCache(nothing, nothing, nothing)

"""
    ThetaLState(φ, L, x0)

State of a curve: φ = θ − α at the N grid points, length L, and x₀ = x(α = 0); `length(X)` is N. A state is a
value: what is computed from it (its spectrum, curvature, points) is kept with it, so never change `φ` in place; make
a new state.
"""
struct ThetaLState
    φ::Vector{Float64}
    L::Float64
    x0::Vec2
    cache::_StateCache
    ThetaLState(φ::Vector{Float64}, L::Float64, x0::Vec2, cache::_StateCache) = new(φ, L, x0, cache)
end
ThetaLState(φ::AbstractVector{<:Real}, L::Real, x0) =
    ThetaLState(convert(Vector{Float64}, φ), Float64(L), convert(Vec2, x0), _StateCache())
Base.length(X::ThetaLState) = length(X.φ)    # N, the number of grid points
Base.show(io::IO, X::ThetaLState) = print(io, "ThetaLState(", X.φ, ", ", X.L, ", ", X.x0, ")")

# the cached quantities (read only: never change what these return)
function _phihat(X::ThetaLState)                 # the half spectrum of φ (unnormalized, as rfft)
    c = X.cache
    c.φh === nothing && (c.φh = _rfft(X.φ))
    c.φh::Vector{ComplexF64}
end
function _theta_alpha(X::ThetaLState)            # θ_α = 1 + φ_α at the grid points
    c = X.cache
    if c.θα === nothing
        N = length(X.φ); w = _work(N); h = _dα_half!(copyto!(_take_h(w), _phihat(X)), N)
        θα = _irfft!(w, Vector{Float64}(undef, N), h); _give!(w, h)
        for j in 1:N
            θα[j] += 1
        end
        c.θα = θα
    end
    c.θα::Vector{Float64}
end
function _points(X::ThetaLState)                 # the grid points
    c = X.cache
    c.points === nothing && (c.points = _curve_points!(Vector{Vec2}(undef, length(X.φ)), X))
    c.points::Vector{Vec2}
end

_α(j::Int, N::Int) = 2π * (j - 1) / N                    # the j-th grid point

"""
Trigonometric interpolant of the curve: returns z(α) as a function (Vec2), using the
closure-corrected Fourier series of x_α = (L/2π)(cos θ, sin θ).
"""
function curve_function(X::ThetaLState)
    N = length(X.φ); sα = X.L / 2π
    zα = Vector{ComplexF64}(undef, N)
    for j in 1:N
        zα[j] = sα * cis(_α(j, N) + X.φ[j])
    end
    c = _fft(zα)
    for j in 1:N
        k = _wavenumber(j, N)                    # the mean (closure) and the Nyquist mode removed
        c[j] = k == 0 || _nyquist(k, N) ? 0.0im : c[j] / (N * im * k)
    end
    _CurveSeries(c, complex(X.x0[1], X.x0[2]) - sum(c), Float64.(_kvec(N)))
end
struct _CurveSeries
    c::Vector{ComplexF64}
    z0::ComplexF64
    k::Vector{Float64}
end
# (callable structs in the short form, forwarding to a function: the long form `function (p::T)(…)` costs a
# second of parser compilation in a fresh Julia 1.12 session, paid by every precompilation)
(z::_CurveSeries)(a::Real) = _series_at(z, Float64(a))
function _series_at(z::_CurveSeries, a::Float64)
    v = z.z0
    @inbounds for i in 2:length(z.c)
        z.c[i] == 0 && continue
        v += z.c[i] * cis(z.k[i] * a)
    end
    Vec2(real(v), imag(v))
end

"""Grid points x(α_j)."""
curve_points(X::ThetaLState) = copy(_points(X))

# the grid points into x: the series of curve_function, summed at the grid points by one inverse FFT (O(N log N))
function _curve_points!(x::Vector{Vec2}, X::ThetaLState)
    N = length(X.φ); w = _work(N); c = _take_c(w); sα = X.L / 2π
    for j in 1:N
        c[j] = sα * cis(_α(j, N) + X.φ[j])
    end
    _fwd!(w, c); z0 = complex(X.x0[1], X.x0[2])
    for j in 1:N
        k = _wavenumber(j, N)                    # the mean (closure) and the Nyquist mode removed
        c[j] = k == 0 || _nyquist(k, N) ? 0.0im : c[j] / (N * im * k)
        z0 -= c[j]
    end
    _bwd!(w, c)
    for j in 1:N
        v = z0 + c[j]
        x[j] = Vec2(real(v), imag(v))
    end
    _give!(w, c)
    x
end

"""Tangent angles θ = α + φ at the grid points."""
function tangent_angles(X::ThetaLState)
    N = length(X.φ); θ = Vector{Float64}(undef, N)
    for j in 1:N
        θ[j] = _α(j, N) + X.φ[j]
    end
    θ
end

"""Curvature κ = θ_s = (1 + φ_α)/s_α at the grid points (spectral; positive where the curve is convex)."""
function curvature(X::ThetaLState)
    θα = _theta_alpha(X); sα = X.L / 2π; κ = Vector{Float64}(undef, length(θα))
    for j in eachindex(κ)
        κ[j] = θα[j] / sα
    end
    κ
end

"""
    arclength_derivative(X, f)

The derivative ∂f/∂s = (2π/L) f_α of a periodic function `f` given at the grid points of `X` (spectral, Nyquist
mode dropped): e.g. κ_ss = `arclength_derivative(X, arclength_derivative(X, curvature(X)))`.
"""
function arclength_derivative(X::ThetaLState, f::AbstractVector{<:Real})
    g = _dα(f); c = 2π / X.L
    for j in eachindex(g)
        g[j] *= c
    end
    g
end

"""
    velocity_components(X, u) -> (U, W)

Normal and tangential parts of the velocities `u` at the grid points of `X`: U = u·n, W = u·τ.
"""
function velocity_components(X::ThetaLState, u::AbstractVector)
    N = length(X.φ); U = Vector{Float64}(undef, N); W = Vector{Float64}(undef, N)
    for j in 1:N
        sn, cs = sincos(_α(j, N) + X.φ[j]); x, y = Float64(u[j][1]), Float64(u[j][2])
        U[j] = sn * x - cs * y; W[j] = cs * x + sn * y
    end
    U, W
end

"""The velocities U n + W τ at the grid points of `X` (the inverse of `velocity_components`)."""
function velocity_vectors(X::ThetaLState, U::AbstractVector, W::AbstractVector)
    N = length(X.φ); u = Vector{Vec2}(undef, N)
    for j in 1:N
        sn, cs = sincos(_α(j, N) + X.φ[j])
        u[j] = Vec2(U[j] * sn + W[j] * cs, -U[j] * cs + W[j] * sn)
    end
    u
end

"""
    curve_moments(X) -> (A, M)

Area A and first moment M = ∫ x dS of the region inside the curve, A = ½∫ x × x_α dα and M = ⅓∫ x (x × x_α) dα by
the trapezoidal rule (spectral).
"""
function curve_moments(X::ThetaLState)
    N = length(X.φ); z = _points(X); sα = X.L / 2π
    A = 0.0; Mx = 0.0; My = 0.0
    for j in 1:N
        sn, cs = sincos(_α(j, N) + X.φ[j])
        w = sα * (z[j][1] * sn - z[j][2] * cs)
        A += w; Mx += z[j][1] * w; My += z[j][2] * w
    end
    A * π / N, Vec2(Mx, My) * (2π / N) / 3
end

"""Area of the region inside the curve."""
area(X::ThetaLState) = curve_moments(X)[1]

"""Centroid of the region inside the curve."""
centroid(X::ThetaLState) = (m = curve_moments(X); m[2] / m[1])

"""Length of the curve."""
perimeter(X::ThetaLState) = X.L

"""Tail of the φ spectrum, the largest |φ̂_k|/N in the band N/2 ± N/8: how well N resolves the curve."""
spectral_tail(X::ThetaLState) = _half_tail(_phihat(X), length(X.φ)) / length(X.φ)

# the φ of a state on twice or half the points (`double_grid`, `halve_grid`, src/grid/fourier.jl)
double_grid(X::ThetaLState) = ThetaLState(_double_grid(_phihat(X), length(X.φ)), X.L, X.x0)
halve_grid(X::ThetaLState) = ThetaLState(_halve_grid(_phihat(X), length(X.φ)), X.L, X.x0)
