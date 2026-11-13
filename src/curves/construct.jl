# States from a curve: from points already at equal arclength (`from_arclength_samples`), and from any
# parametrization by Fourier series (`from_curve`).

# φ = θ − α from tangent angles θ at the N grid points: θ made continuous first
function _φ_from_angles!(θ::Vector{Float64})
    N = length(θ)
    for j in 2:N
        θ[j] += 2π * round((θ[j-1] - θ[j]) / 2π)
    end
    for j in 1:N
        θ[j] -= _α(j, N)
    end
    θ
end

"""
    from_arclength_samples(points) -> ThetaLState

The state of the closed counterclockwise curve through `points`: N points equispaced in arclength, the first at
α = 0. The tangent comes from the spectral derivative of the trigonometric interpolant of the points.
"""
function from_arclength_samples(points::AbstractVector)
    N = length(points); z = Vector{ComplexF64}(undef, N)
    for j in 1:N
        z[j] = complex(Float64(points[j][1]), Float64(points[j][2]))
    end
    zh = _fft(z)
    for j in 1:N
        k = _wavenumber(j, N)
        zh[j] = _nyquist(k, N) ? 0.0im : im * k * zh[j]
    end
    zα = _ifft(zh); θ = Vector{Float64}(undef, N); L = 0.0
    for j in 1:N
        θ[j] = angle(zα[j]); L += abs(zα[j])
    end
    ThetaLState(_φ_from_angles!(θ), 2π * L / N, Vec2(real(z[1]), imag(z[1])))
end

# Fourier data of a parametrized curve z(t) = γ(t): z_t and the arclength s(t) = (L/2π) t + P(t) − P(0), P the
# periodic antiderivative of the speed − L/2π, with the arclength sg on the M grid.
struct _ArcMap
    k::Vector{Float64}
    zth::Vector{ComplexF64}
    sh::Vector{ComplexF64}
    P0::Float64
    L::Float64
    sg::Vector{Float64}
end

# M is doubled from 64 until the tails of ẑ and of the speed are below tol/100, or stop falling below 1e-10
# (round-off: a smaller target sends M to the cap, and sums over that many terms lose 1e-10 in φ).
Base.@nospecializeinfer function _arclength_map(γ, tol::Float64; M::Int=64, Mmax::Int=2^18)
    @nospecialize γ                                         # one compiled map for every curve
    prev = Inf
    while true
        z = Vector{ComplexF64}(undef, M)
        for j in 1:M
            x = γ(_α(j, M)); z[j] = complex(Float64(x[1]), Float64(x[2]))
        end
        zh = _fft(z); zth = similar(zh)
        for j in 1:M
            k = _wavenumber(j, M); zh[j] /= M
            zth[j] = _nyquist(k, M) ? 0.0im : im * k * zh[j]
        end
        zt = _ifft(zth); v = Vector{Float64}(undef, M)
        for j in 1:M
            v[j] = abs(zt[j]) * M                           # the speed |z_t|
        end
        vh = _fft(v)
        for j in 1:M
            vh[j] /= M
        end
        tail = max(_rel_tail(zh), _rel_tail(vh))
        if tail < 1e-2 * tol || (tail < 1e-10 && tail > prev / 4) || 2M > Mmax
            L = 2π * real(vh[1]); sh = similar(vh)
            for j in 1:M
                k = _wavenumber(j, M)
                sh[j] = k == 0 || _nyquist(k, M) ? 0.0im : vh[j] / (im * k)
            end
            P = _ifft(sh); P0 = real(P[1]) * M; sg = Vector{Float64}(undef, M)
            for j in 1:M
                sg[j] = L * (j - 1) / M + real(P[j]) * M - P0
            end
            return _ArcMap(Float64.(_kvec(M)), zth, sh, P0, L, sg)
        end
        M *= 2; prev = tail
    end
end

# s(t) and z_t(t) from the series
function _s_zt(m::_ArcMap, t::Float64)
    s = 0.0; zt = 0.0im
    @inbounds for i in eachindex(m.k)
        e = cis(m.k[i] * t)
        s += real(m.sh[i] * e); zt += m.zth[i] * e
    end
    m.L * t / 2π + s - m.P0, zt
end

# the point γ(0), and the state with N points at equal arclength: Newton on s(t) = jL/N from the linear inverse of s on
# the grid
_start(@nospecialize(γ)) = (x = γ(0.0); Vec2(x[1], x[2]))
function _from_map(m::_ArcMap, x0::Vec2, N::Int)
    θ = zeros(N)
    Threads.@threads for j in 1:N
        θ[j] = _angle_at(m, (j - 1) * m.L / N)
    end
    ThetaLState(_φ_from_angles!(θ), m.L, x0)
end
function _angle_at(m::_ArcMap, sj::Float64)
    M = length(m.k)
    i = clamp(searchsortedlast(m.sg, sj), 1, M)
    t = 2π * (i - 1) / M + (sj - m.sg[i]) / ((i == M ? m.L : m.sg[i+1]) - m.sg[i]) * (2π / M)
    for _ in 1:30
        s, zt = _s_zt(m, t)
        δ = (s - sj) / abs(zt); t -= δ
        abs(δ) < 1e-14 && break
    end
    angle(_s_zt(m, t)[2])
end

"""
    from_curve(γ, N; tol=1e-6) -> ThetaLState
    from_curve(γ; tol=1e-6, N=64, Nmax=16384) -> ThetaLState

The state of the closed counterclockwise curve ξ ↦ γ(ξ) (2π-periodic; a point is anything indexed by 1, 2) with N
points at equal arclength, the first at γ(0). Spectral: γ is sampled at M equispaced ξ, M doubled until the
Fourier tails of γ and of its speed are below tol/100; the arclength is the antiderivative of the speed, the
points are found by Newton on it, and θ is the angle of γ_ξ there. Without N, N is doubled from `N` until the tail
of the φ spectrum (`spectral_tail`) is below tol/100, or until doubling would exceed `Nmax`.
"""
from_curve(@nospecialize(γ), N::Int; tol::Real=1e-6) =
    _from_map(_arclength_map(γ, Float64(tol)), _start(γ), N)
from_curve(@nospecialize(γ); tol::Real=1e-6, N::Int=64, Nmax::Int=16384) = _from_curve(γ, Float64(tol), N, Nmax)
Base.@nospecializeinfer function _from_curve(@nospecialize(γ), tol::Float64, N::Int, Nmax::Int)
    m = _arclength_map(γ, tol); x0 = _start(γ)
    while true
        X = _from_map(m, x0, N)
        (spectral_tail(X) < 1e-2 * tol || 2N > Nmax) && return X
        N *= 2
    end
end
