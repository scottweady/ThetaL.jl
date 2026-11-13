# θ–L kinematics (Hou, Lowengrub & Shelley 1994): the time derivatives of a state from its normal velocity, with
# the tangential velocity that keeps the points at equal arclength (equations and conventions: the module
# docstring); the rigid motion in a normal velocity (`rigid_part`), which the step takes out and integrates
# exactly; and the velocity of the tracked points (`frame_velocity`).

# The tangential velocity T of Hou, Lowengrub & Shelley for the normal velocity Uf (already filtered): equal
# arclength kept, T(0) = 0, written into T (one real FFT and its inverse); returns L_t. θ_α from the state's cache.
function _hls_tangential_velocity!(T::Vector{Float64}, X::ThetaLState, Uf::Vector{Float64})
    N = length(X.φ); w = _work(N); θα = _theta_alpha(X)
    Lt = 0.0
    for j in 1:N
        Lt += Uf[j] * θα[j]
    end
    Lt *= 2π / N
    f = _take_r(w)
    for j in 1:N
        f[j] = Lt / 2π - Uf[j] * θα[j]
    end
    _intα!(T, f); _give!(w, f)
    Lt
end

"""
    thetaL_rhs(X, U, filt) -> (θ_t, L_t, x₀_t)

Time derivatives of the state `X` under the normal velocity `U` (filtered by the mode weights `filt` first): θ_t on
the grid, L_t, and the velocity of the point x₀. No rigid part is taken out here (the step does that).
"""
function thetaL_rhs(X::ThetaLState, U::AbstractVector, filt::AbstractVector)
    N = length(X.φ); θth = Vector{ComplexF64}(undef, N ÷ 2 + 1)
    Lt, drift = _thetaL_rhs_hat!(θth, X, _vec64(U), Vector{Float64}(filt[1:N÷2+1]))   # filt on the half spectrum
    _irfft!(_work(N), Vector{Float64}(undef, N), θth), Lt, drift
end

# The half spectrum of θ_t = (T θ_α − U_α)/s_α into θth (and T into Td, if given); returns L_t and the velocity of
# x₀. Five real FFTs: U, the filtered U back on the grid, the tangential velocity T (two), and T θ_α; U_α only as
# the spectrum ik·filt·Û, and θ_α from the state's cache.
function _thetaL_rhs_hat!(θth::Vector{ComplexF64}, X::ThetaLState, U::Vector{Float64}, filt::Vector{Float64},
                          Td::Union{Nothing,Vector{Float64}}=nothing)
    N = length(X.φ); sα = X.L / 2π; w = _work(N)
    Uh = _rfft!(w, _take_h(w), U); c = _take_h(w)
    for j in eachindex(Uh)
        c[j] = filt[j] * Uh[j]
    end
    Uf = _irfft!(w, _take_r(w), c)
    T = _take_r(w); Lt = _hls_tangential_velocity!(T, X, Uf); θα = _theta_alpha(X)
    Td === nothing || copyto!(Td, T)                     # the deformation's tangential velocity, if asked for
    for j in 1:N
        T[j] *= θα[j]
    end
    _rfft!(w, θth, T)
    for j in eachindex(θth)
        k = j - 1
        Uαh = _nyquist(k, N) ? 0.0im : im * k * filt[j] * Uh[j]
        θth[j] = (θth[j] - Uαh) / sα
    end
    sn, cs = sincos(X.φ[1]); drift = Uf[1] * Vec2(sn, -cs)
    _give!(w, Uh, c, Uf, T)
    Lt, drift
end

"""
    RigidMotion(translation, rotation_rate, center)

A rigid velocity field u(x) = translation + rotation_rate · ẑ × (x − center) (`rigid_velocity(m, x)`).
"""
struct RigidMotion
    translation::Vec2
    rotation_rate::Float64
    center::Vec2
end
const NO_RIGID_MOTION = RigidMotion(Vec2(0.0, 0.0), 0.0, Vec2(0.0, 0.0))

"""The velocity of the rigid motion `m` at the point `x`."""
function rigid_velocity(m::RigidMotion, x)
    r = Vec2(x[1], x[2]) - m.center
    m.translation + m.rotation_rate * Vec2(-r[2], r[1])
end

"""
    rigid_part(X, U; mode = :full) -> (motion::RigidMotion, U_deformation)

The rigid motion in the normal velocity `U` of the curve `X`, and the rest. The motion is the translation V and
rotation rate Ω about the center x_c (the mean of the grid points) whose normal velocity
U_rig = V·n − Ω (x − x_c)·τ is closest to `U` on the grid (least squares; equal arclength);
`U_deformation = U − U_rig`.

- `mode = :full`: translation and rotation. On a circle a rotation has no normal velocity, so on a near-circle the
  fit turns small errors in `U` into a large Ω: up to 1/(r √ρ) per unit rms of `U`, with r the rms radius,
  ρ = Σ p² / Σ |x − x_c|² over the grid and p the part of (x − x_c)·τ that no translation gives (ρ ≈ ε² for
  an ellipse of axes 1 and 1 + ε, ≈ 0.2 for a two-lobed shape). So the rotation is weighted by a C¹ ramp in ρ:
  the plain least-squares Ω for ρ ≥ 10⁻³ (rigid rotation exact), none for ρ ≤ 10⁻⁴, the smoothstep between; the
  gain stays below 33/r. The rotation not taken out stays in `U_deformation`, which the step integrates (correct,
  not exact at any Δt; a near-circle's rotation hardly changes it).
- `mode = :translation`: translation only.
- `mode = :none`: no rigid part (`NO_RIGID_MOTION`, U returned as it is).
"""
function rigid_part(X::ThetaLState, U::AbstractVector; mode::Symbol=:full)
    Ud = Vector{Float64}(undef, length(X.φ))
    _rigid_part!(Ud, X, U, mode), Ud
end

# The weight of the fitted rotation for ρ = Σp²/Σ|x − x_c|² (`rigid_part`): 0 for ρ ≤ ρ₀, 1 for ρ ≥ ρ₁, the C¹
# smoothstep between (Ω continuous and C¹ as a shape crosses the band; a hard cutoff would jump by 1/(r√ρ) times
# the error in U). The noise gain s(ρ)/(r√ρ) is at most 32.2/r. Fixed: ρ is a shape measure (free of scale and N),
# and only the exactness of a near-circle's rotation, which barely moves it, depends on the band.
const _ROTATION_RAMP = (1e-4, 1e-3)
function _rotation_weight(ρ::Float64)
    ρ₀, ρ₁ = _ROTATION_RAMP
    ρ >= ρ₁ && return 1.0
    ρ <= ρ₀ && return 0.0
    t = (ρ - ρ₀) / (ρ₁ - ρ₀)
    t^2 * (3 - 2t)
end

# the rigid motion, with U_deformation written into Ud
function _rigid_part!(Ud::Vector{Float64}, X::ThetaLState, U::AbstractVector, mode::Symbol)
    mode in (:full, :translation, :none) || throw(ArgumentError("rigid mode must be :full, :translation or :none, got :$mode"))
    N = length(X.φ)
    if mode === :none
        for j in 1:N
            Ud[j] = U[j]
        end
        return NO_RIGID_MOTION
    end
    w = _work(N); x = _points(X)
    xc = Vec2(0.0, 0.0)
    for j in 1:N
        xc += x[j]
    end
    xc /= N
    sn = _take_r(w); cs = _take_r(w); b = _take_r(w)                    # n = (sin θ, −cos θ); b = −(x − xc)·τ
    r2 = 0.0
    for j in 1:N
        sn[j], cs[j] = sincos(_α(j, N) + X.φ[j]); r = x[j] - xc
        b[j] = -(r[1] * cs[j] + r[2] * sn[j]); r2 += r[1]^2 + r[2]^2
    end
    # translation: the 2×2 normal equations for n_x = sin θ, n_y = −cos θ, applied to U and to b
    g11 = 0.0; g12 = 0.0; g22 = 0.0; u1 = 0.0; u2 = 0.0; b1 = 0.0; b2 = 0.0
    for j in 1:N
        g11 += sn[j]^2; g12 -= sn[j] * cs[j]; g22 += cs[j]^2
        u1 += sn[j] * U[j]; u2 -= cs[j] * U[j]; b1 += sn[j] * b[j]; b2 -= cs[j] * b[j]
    end
    det = g11 * g22 - g12^2                         # (no closure here: it would box the sums above)
    tU = Vec2((g22 * u1 - g12 * u2) / det, (g11 * u2 - g12 * u1) / det)
    tb = Vec2((g22 * b1 - g12 * b2) / det, (g11 * b2 - g12 * b1) / det)
    # rotation: U against the part of b that no translation gives
    Ω = 0.0
    if mode === :full
        bb = 0.0; bU = 0.0
        for j in 1:N
            p = b[j] - (tb[1] * sn[j] - tb[2] * cs[j])
            bb += p^2; bU += p * U[j]
        end
        s = _rotation_weight(bb / r2)                  # both sums over the grid: ρ independent of N
        Ω = s > 0 ? s * bU / bb : 0.0                  # (the rest of the rotation stays in Ud, by V below)
    end
    V = tU - Ω * tb
    for j in 1:N
        Ud[j] = U[j] - (V[1] * sn[j] - V[2] * cs[j] + Ω * b[j])
    end
    _give!(w, sn, cs, b)
    RigidMotion(V, Ω, xc)
end

"""
    frame_velocity(X, U; rigid = :full) -> Vector{Vec2}

The velocity of the tracked grid points of `X` under the normal velocity `U`: U n + T τ, with T the tangential
velocity of the step (the rigid motion's, plus the equal-arclength one of the deformation, T(0) = 0). This is what
a moving bulk mesh attached to the curve's points follows. `rigid` as in `rigid_part`.
"""
function frame_velocity(X::ThetaLState, U::AbstractVector; rigid::Symbol=:full)
    N = length(X.φ)
    motion, Ud = rigid_part(X, U; mode=rigid)
    Td = Vector{Float64}(undef, N); _hls_tangential_velocity!(Td, X, Ud)
    x = curve_points(X); u = Vector{Vec2}(undef, N)
    for j in 1:N
        sn, cs = sincos(_α(j, N) + X.φ[j])
        n = Vec2(sn, -cs); τ = Vec2(cs, sn)
        u[j] = rigid_velocity(motion, x[j]) + Ud[j] * n + Td[j] * τ
    end
    u
end
