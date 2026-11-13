# The stiffness of a velocity, measured (`measure_symbol`), and the symbol given for it checked against it
# (`check_symbol`). For a velocity whose small-scale symbol is not known in closed form (a boundary integral, a
# numerical solve), or to confirm one that is.
#
# Measurement: the linearization J of φ ↦ θ_t (the step's kinematics: rigid part out, Hou–Lowengrub–Shelley θ_t of
# the deformation, rotation rate added to every angle) at the state X of one curve, the other curves held, applied to
# packets p = w(α) e^{ik(α−α₀)}: J p = J(w cos) + i J(w sin) by linearity, each from a central difference
# (θ_t(φ + εu) − θ_t(φ − εu))/2ε, which cancels the even nonlinear terms (O(ε²) error). The rate is the Rayleigh
# quotient λ = ⟨p, J p⟩/⟨p, p⟩ on the grid: with w ≡ 1 (one window) the Fourier coefficient (J e^{ikα})_k/N, exact for a
# circle and the mode-averaged rate on other curves; with m windows, w a periodic Gaussian about α₀ = 2π(j−1)/m, the
# local rate there (the stiffness of U = −a(α)κ is a(α₀)q², which the global mode averages). Four evaluations per k
# and window.
#
# Stability: the step treats σ_g (the given symbol) exactly and R = (λ − σ_g)φ̂ explicitly. Per mode the scheme maps
# φ̂ ↦ G φ̂, G from its weights (`scheme_weights`) at z = σ_g h, with the explicit rate μ = λ − σ_g. For ETDRK2,
#   G = e^{z} + hμ(ϕ₁ − 2ϕ₂)(z) + 2hμ ϕ₂(z) (e^{z/2} + (hμ/2) ϕ₁(z/2)),
# and as σ_g h → −∞, G → r + 2r² with r = μ/σ_g = λ/σ_g − 1: |G| ≤ 1 for every Δt only if −1 ≤ r ≤ 1/2, i.e.
# 0 ≤ λ/σ_g ≤ 3/2 for real rates (ExpEuler: G → −r, 0 ≤ λ/σ_g ≤ 2). A symbol weaker than 2/3 of the true rate (or
# none) leaves a step limit Δt ~ 1/|λ| at that mode. The check computes, per measured mode, the largest step with
# |G| ≤ max(1, e^{Re λ h}) (a mode that truly grows may grow as fast as it does), by scanning h and bisecting.

"""
    measure_symbol(problem, states = nothing; curve = 1, k = nothing, ε = 1e-4, windows = 1, t = 0.0,
                   warm = nothing, tol = 1e-6, Nmax = 16384)
    measure_symbol(integ; curve = 1, k = nothing, ε = 1e-4, windows = 1)
        -> (; k, q, σ, λ, given)

The growth rate of small-scale modes of curve `curve` under the problem's velocity, measured: φ of that curve is
perturbed by ±ε w(α) cos(k(α − α₀)) and ±ε w(α) sin(k(α − α₀)) (L, x₀ and the other curves held), θ_t is computed
from the velocity as the step does (the rigid part taken out per `problem.rigid`, no filter), and the rate
λ = ⟨p, J p⟩/⟨p, p⟩ of the packet p = w e^{ik(α−α₀)} is read from the central differences (four evaluations per k
and window). `states` default to the problem's curves (parametrizations resolved to `tol`, as `init` does); the
integrator's method measures at its time, states, context and carried state.

- `k`: the Fourier modes (default: six from 2 to N/4, about geometrically spaced);
- `ε`: the amplitude of the perturbation of φ (radians). The central difference has an O(ε²) error; noise η in U
  enters λ as about q η/ε: raise ε for a velocity computed to a loose tolerance;
- `windows`: `1`, the whole curve (w ≡ 1: λ is the Fourier coefficient (J e^{ikα})_k/N, the mode-averaged rate);
  `m > 1`, m Gaussian packets centered at α₀ = 2π(j − 1)/m, of width π/(m√2) in α: the local rate, which is what a
  frozen symbol must bound on a curve whose stiffness varies along it (anisotropy);
- `t`, `warm`: the evaluations' time and the warm start of the first (each passes its `ctx` to the next).

Returns `k`, `q = 2π k/L`, the measured `σ` (Re λ, the most negative over the windows), `λ` (complex, one column
per window: an imaginary part is a transport rate, which the real symbol leaves explicit), and `given` (the
problem's symbol at q, zero without one). The velocity is told `ev.kind = :probe`.
"""
Base.@nospecializeinfer function measure_symbol(problem::CurveProblem, states=nothing; curve::Int=1, k=nothing,
                                                ε::Real=1e-4, windows::Int=1, t::Real=0.0, warm=nothing,
                                                tol::Real=1e-6, Nmax::Int=16384)
    @nospecialize states warm k
    Xs = states === nothing ?
        ThetaLState[c isa ThetaLState ? c : _from_curve(c, Float64(tol), 64, Nmax) for c in problem.curves] :
        Vector{ThetaLState}(collect(states))
    n = length(Xs)
    fields = isempty(problem.fields) ? _NO_FIELDS :
        Base.inferencebarrier(_initial_fields_tuple)(problem, Xs)::NamedTuple
    bulk = problem.bulk === nothing ? nothing : copy(problem.bulk.y0)
    _measure_symbol(problem, Xs, curve, k, Float64(ε), windows, Float64(t), warm, collect(1:n),
                    Any[c <= length(problem.tags) ? problem.tags[c] : nothing for c in 1:n], 0, fields, bulk)
end

function measure_symbol(integ::ThetaLIntegrator; curve::Int=1, k=nothing, ε::Real=1e-4, windows::Int=1)
    @nospecialize k
    _measure_symbol(integ.problem, integ.states, curve, k, Float64(ε), windows, integ.t, integ.ctx, integ.ids,
                    integ.tags, integ.topology, integ.fields, integ.bulk)
end

# the default modes: about geometric from 2 to N/4
_probe_modes(N::Int) = unique(round.(Int, exp.(range(log(2), log(max(N ÷ 4, 2)), length=6))))

Base.@nospecializeinfer function _measure_symbol(problem::CurveProblem, states::Vector{ThetaLState}, curve::Int, k,
                                                 ε::Float64, windows::Int, t::Float64, warm, ids::Vector{Int},
                                                 tags::Vector{Any}, topology::Int, fields::NamedTuple,
                                                 bulk::Union{Nothing,Vector{Float64}})
    @nospecialize k warm
    1 <= curve <= length(states) || throw(ArgumentError("curve $curve: the problem has $(length(states)) curves"))
    windows >= 1 || throw(ArgumentError("windows must be at least 1"))
    ε > 0 || throw(ArgumentError("ε must be positive"))
    X = states[curve]; N = length(X.φ)
    ks = k === nothing ? _probe_modes(N) : collect(Int, k)
    all(kk -> 1 <= kk < N ÷ 2, ks) || throw(ArgumentError("the modes k must be in 1:$(N ÷ 2 - 1)"))
    α = [_α(j, N) for j in 1:N]; ones_f = ones(N)
    ev = Evaluation(t, :probe, 0, 0.0, 0.0, 1, ids, tags, topology, states, fields, bulk)
    ctx = Ref{Any}(warm)
    # θ_t of the curve with φ + δ, as the step forms it from the velocity
    function θt(δ::Vector{Float64})
        Xp = ThetaLState(X.φ .+ δ, X.L, X.x0)
        sp = copy(states); sp[curve] = Xp
        U, c, _ = _evaluate_velocity(problem, sp, ev, ctx[]); ctx[] = c
        motion, Ud = rigid_part(Xp, U[curve]; mode=problem.rigid)
        thetaL_rhs(Xp, Ud, ones_f)[1] .+ motion.rotation_rate
    end
    Jd(u) = (θt(ε .* u) .- θt(-ε .* u)) ./ (2ε)          # J u by a central difference
    λ = Matrix{ComplexF64}(undef, length(ks), windows)
    for (i, kk) in enumerate(ks), j in 1:windows
        α0 = 2π * (j - 1) / windows
        w = windows == 1 ? ones_f : [exp(2 * (windows / π)^2 * (cos(a - α0) - 1)) for a in α]
        uc = w .* cos.(kk .* (α .- α0)); us = w .* sin.(kk .* (α .- α0))
        Jp = Jd(uc) .+ im .* Jd(us)                       # J p, p = uc + i us = w e^{ik(α−α₀)}
        λ[i, j] = sum(conj.(uc .+ im .* us) .* Jp) / sum(abs2, w)
    end
    q = (2π / X.L) .* ks
    (; k=ks, q, σ=vec(minimum(real.(λ); dims=2)), λ, given=_curve_rates(problem.symbol, X, ks))
end

"""
    check_symbol(problem, states = nothing; dt = nothing, scheme = ETDRK2(), curves = nothing, warn = true, kw...)
    check_symbol(integ; dt = nothing, curves = nothing, warn = true, kw...)
        -> (; ok, dt_limit, curves)

Whether the problem's symbol bounds the stiffness of its velocity: the rates are measured (`measure_symbol`; `kw`
are passed on: `k`, `ε`, `windows`, …) on each curve of `curves` (default all), and for each mode the largest step
`dt_max` that the scheme takes stably is computed, with the symbol σ_g treated exactly and the rest,
(λ − σ_g)φ̂, explicitly: per mode the step multiplies φ̂ by G(σ_g h, (λ − σ_g)h) (from the scheme's weights), and
`dt_max` is the first h at which |G| exceeds max(1, e^{Re λ h}) (`Inf`: stable at every step).

For ETDRK2 and real rates the condition at large steps is 0 ≤ λ/σ_g ≤ 3/2: a symbol may overstate the true rate
(the mode is then damped less per step), and must not understate it by more than a third (ExpEuler:
0 ≤ λ/σ_g ≤ 2, by up to half). Beyond that a step limit
Δt ~ 1/|λ| remains at that mode, which a stiff velocity turns into Δt ~ q^{−p}.

`ok`: with `dt`, whether `dt ≤ dt_limit`; without, whether no measured mode has a limit (`dt_limit = Inf`). A
finite limit at low q only can be the explicit treatment of a part that is not stiff (pass the step you plan as
`dt`). `curves[c]` is the measurement of curve c with `dt_max`. With `warn`, a failed check warns with the worst
modes. Only the measured modes are checked: the ratio λ/σ_g at the largest measured q stands for the higher ones.
"""
function check_symbol(problem::CurveProblem, states=nothing; curves=nothing, kw...)
    @nospecialize states curves
    n = states === nothing ? length(problem.curves) : length(states)
    m = Base.structdiff(values(kw), NamedTuple{(:dt, :scheme, :warn)})
    _check(c -> measure_symbol(problem, states; curve=c, m...), curves === nothing ? (1:n) : curves; kw...)
end
function check_symbol(integ::ThetaLIntegrator; curves=nothing, kw...)
    @nospecialize curves
    m = Base.structdiff(values(kw), NamedTuple{(:dt, :scheme, :warn)})
    _check(c -> measure_symbol(integ; curve=c, m...), curves === nothing ? eachindex(integ.states) : curves;
           scheme=integ.scheme, kw...)
end

function _check(measure, curves; dt=nothing, scheme::OneStepScheme=ETDRK2(), warn::Bool=true, kw...)
    res = Any[]; lim = Inf
    for c in curves
        r = measure(c)
        dmax = [minimum(_stable_step(scheme, r.given[i], r.λ[i, j]) for j in axes(r.λ, 2)) for i in eachindex(r.k)]
        push!(res, merge(r, (; curve=c, dt_max=dmax))); lim = min(lim, minimum(dmax))
    end
    ok = dt === nothing ? isinf(lim) : dt <= lim
    if warn && !ok
        lines = String[]
        for r in res, i in eachindex(r.k)
            isfinite(r.dt_max[i]) && (dt === nothing || r.dt_max[i] < dt) &&
                push!(lines, "curve $(r.curve), q = $(round(r.q[i]; sigdigits=4)): measured σ = " *
                             "$(round(r.σ[i]; sigdigits=4)), symbol $(round(r.given[i]; sigdigits=4)), " *
                             "stable for Δt < $(round(r.dt_max[i]; sigdigits=3))")
        end
        @warn "the symbol does not bound the measured stiffness" *
              (dt === nothing ? "" : " at Δt = $dt") * ":\n  " * join(lines, "\n  ")
    end
    (; ok, dt_limit=lim, curves=res)
end

# the amplification of one mode over a step h: the scheme with z = σh exact and the rate μ explicit,
# y' = σy + μy (stages y_i = e^{c_i z} y + h Σ_j A_ij μ y_j, end e^{z} y + h Σ_i b_i μ y_i)
function _amplification(scheme::OneStepScheme, σ::Float64, μ::ComplexF64, h::Float64)
    z = σ * h; nodes = stage_nodes(scheme); A, b, _ = scheme_weights(scheme, z)
    y = Vector{ComplexF64}(undef, length(nodes))
    for i in eachindex(nodes)
        s = complex(exp(nodes[i] * z))
        for j in 1:i-1
            s += h * A[i][j] * μ * y[j]
        end
        y[i] = s
    end
    g = complex(exp(z))
    for i in eachindex(nodes)
        g += h * b[i] * μ * y[i]
    end
    g
end

# the largest step at which the mode stays bounded, |G| ≤ max(1, e^{Re λ h}): scanned over h from 10⁻³ to 10⁸
# times 1/max(|σ|, |λ|), the first failure refined by bisection; Inf if none
function _stable_step(scheme::OneStepScheme, σ::Float64, λ::ComplexF64)
    s = max(abs(σ), abs(λ)); s > 0 || return Inf
    μ = λ - σ
    bad(h) = abs(_amplification(scheme, σ, μ, h)) > max(1.0, exp(real(λ) * h)) * (1 + 1e-10)
    hs = exp10.(range(-3, 8, length=221)) ./ s
    i = findfirst(bad, hs)
    i === nothing && return Inf
    i == 1 && return hs[1]
    lo, hi = hs[i-1], hs[i]
    for _ in 1:50
        mid = sqrt(lo * hi)
        bad(mid) ? (hi = mid) : (lo = mid)
    end
    lo
end
