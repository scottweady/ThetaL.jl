# The drivers: run an integrator over a time span with a fixed step (or a step function) or a step chosen from the
# embedded error estimate (`adaptive`), landing on the saved times, adapting N, and keeping the curves at the saved
# times. `evolve` for a curve (or several) and a simple
# velocity U(X, t); `solve` for a `CurveProblem`. Their result (`Evolution`, `sol(t)`) and a CSV writer for it
# (`save_curves`).

"""
    Evolution

The result of `evolve` and `solve`: the saved times `t`, the curves there (`curves`: a `ThetaLState` each for
`evolve` of one curve, otherwise a vector of states per time), the curve `ids` at each saved time, the number of
accepted `steps` and of rejected ones (`rejections`, with `adaptive`); with state carried with the curves, `fields` (per saved time a named tuple, one vector per curve
each: `sol.fields[i].c`) and `bulk` (per saved time). `sol(t)` is the curves saved at time `t`. With DotPlots
loaded, `dotplot(sol)` draws them in the terminal, one frame per saved time.
"""
struct Evolution{C}
    t::Vector{Float64}
    curves::Vector{C}
    ids::Vector{Vector{Int}}
    steps::Int
    rejections::Int
    fields::Vector{NamedTuple}
    bulk::Vector{Union{Nothing,Vector{Float64}}}
end

# (callable structs in the short form, forwarding to a function: the long form `function (p::T)(…)` costs a
# second of parser compilation in a fresh Julia 1.12 session, paid by every precompilation)
(sol::Evolution)(t::Real) = _saved_curves(sol, Float64(t))
function _saved_curves(sol::Evolution, t::Float64)
    for i in eachindex(sol.t)
        abs(sol.t[i] - t) <= 1e-12 * (1 + max(abs(sol.t[i]), abs(t))) && return sol.curves[i]
    end
    throw(ArgumentError("no curves were saved at t = $t"))
end

# `evolve`'s simple velocity U(X, t) (one curve) or U(Xs, t) (several) as the velocity of a CurveProblem
struct _SimpleVelocity
    U::Any
    several::Bool
end
(f::_SimpleVelocity)(states::Vector{ThetaLState}, ev::Evaluation; warm=nothing) = _simple_velocity(f, states, ev)
function _simple_velocity(f::_SimpleVelocity, states::Vector{ThetaLState}, ev::Evaluation)
    r = f.several ? f.U(states, ev.t) : f.U(states[1], ev.t)
    r isa Tuple && throw(ArgumentError("U returned a Tuple: return the normal velocity only (for a velocity field u, " *
                                       "velocity_components(X, u)[1])"))
    f.several ? r : [r]
end

# the times to stop at in (t0, t1], increasing, t1 last
Base.@nospecializeinfer function _stops(saveat, t0::Float64, t1::Float64)
    @nospecialize saveat
    stops = Float64[]
    if saveat isa Real
        for t in t0+saveat:saveat:t1
            push!(stops, Float64(t))
        end
    elseif saveat !== nothing
        for t in saveat
            τ = Float64(t)
            t0 < τ <= t1 && insert!(stops, searchsortedfirst(stops, τ), τ)
        end
    end
    (isempty(stops) || last(stops) != t1) && push!(stops, t1)          # the end is always reached (and kept)
    stops
end

"""
    evolve(curve, U, tspan; dt, symbol = nothing, saveat = nothing, tol = 1e-6, Nmax = 16384, scheme = ETDRK2(),
           rigid = :full, callback = nothing, coarsen = false, verbose = false,
           adaptive = false, errtol = 1e-5, dtmin = 0.0, dtmax = Inf) -> Evolution

Move the closed counterclockwise `curve` over `tspan = (t0, t1)` with outward normal velocity `U`. One call for
simple scripts; `init`/`advance!` for anything more.

- `curve`: a parametrization `ξ -> (x, y)` (2π-periodic), a `ThetaLState`, or a vector of either: several curves
  that move together.
- `U(X, t)`: the normal velocity at the N points of the curve `X` (a `ThetaLState`: use `curve_points(X)`,
  `curvature(X)`, `tangent_angles(X)`, …; for a velocity field u, `velocity_components(X, u)[1]`); for several
  curves `U(Xs, t)`, one vector per curve. Its rigid part is integrated exactly (`rigid`, as in `CurveProblem`).
- `symbol`: the growth rate σ of a small-scale mode, θ_t ≈ σ φ̂, as `σ(q)` of the wavenumber q = 2π|k|/L along
  the curve or `σ(q, X)`: `q -> -q^2` for U = −κ, `q -> -q^4` for U = κ_ss; `nothing`: none.
- `dt`: the step, a number or a function `dt(X, t)` (`dt(Xs, t)` for several curves). With `adaptive`, optional:
  a number is the first step (default: a guess from the velocity), a function a bound on each step (a stability
  limit, as `dtmax`).
- `adaptive = true`: the step is chosen from the embedded error estimate of the scheme (`ETDRK2`; an
  `ArgumentError` for a scheme without one) and kept below `errtol` (`advance!`'s, an absolute bound per step on
  the largest of |Δφ|, |ΔL|/L and |Δx₀|/L between the step and its embedded first-order one, over the curves).
  The next step is h·clamp(0.9 (errtol/err)^½, 0.2, 5) (err ∝ h²); a step over `errtol` is rejected before its end
  evaluation (one evaluation lost) and retried with the step so shrunk, a step rejected for another reason (the
  velocity threw `StepFailure`, a non-positive length, `validate`) with a quarter of it. The steps stay within
  `dtmin` and `dtmax`; below `dtmin`, or after 30 rejections in a row, the run stops with an error. Small steps
  where the curve moves fast (an initial layer), large ones where it relaxes.
- `saveat`: the times to keep: `nothing` (the start and the end), a spacing, or a list. Steps land on them exactly
  (with `adaptive`, the last step before a saved time is stretched by up to 10% or the last two split evenly, so
  that no sliver is left, and the step from before resumes after it).
- `tol`: N is the smallest power of two from 64 whose φ tail (`spectral_tail`) is below tol/100; it is doubled (up
  to `Nmax`) after any step whose tail exceeds `tol`, and with `coarsen = true` halved when the tail on half the
  points has stayed below tol/100 for 10 steps. Without `coarsen`, N never comes back down: one step of noise in a
  computed velocity is enough to double it for the rest of the run (`ETDRK2(filter_cutoff = :auto)` damps that
  noise). Set `Nmax` to what your velocity can afford: for an O(N²) solve the default 16384 is out of reach.
- `callback(integ)`: called after every accepted step; return `:stop` to end the run. It may edit the curves.
- `verbose`: print a line at the start and at each saved time (`progress_line`: t, steps, rejections, the mean Δt
  since the previous line, N, evaluations), with the wall time so far, and a line whenever N changes (when, and
  why).

Without `adaptive` a rejected step is an error (the step is given: there is nothing to retry with); to control
steps otherwise, use the integrator.
"""
function evolve(curve, U, tspan; dt=nothing, symbol=nothing, saveat=nothing, tol::Real=1e-6, Nmax::Int=16384,
                scheme::OneStepScheme=ETDRK2(), rigid::Symbol=:full, callback=nothing, coarsen::Bool=false,
                verbose::Bool=false, adaptive::Bool=false, errtol::Real=1e-5, dtmin::Real=0.0, dtmax::Real=Inf)
    ctrl = _step_control(dt, scheme, adaptive, errtol, dtmin, dtmax)
    _evolve(curve, U, Float64(tspan[1]), Float64(tspan[2]), dt, symbol, saveat, Float64(tol), Nmax, scheme, rigid,
            callback, coarsen, verbose, ctrl)
end
# (the caller's functions positional and unspecialized: one compiled `evolve` for every script, whose own closures
# cost only the keyword wrapper above)
Base.@nospecializeinfer function _evolve(curve, U, t0::Float64, t1::Float64, dt, symbol, saveat, tol::Float64,
                                         Nmax::Int, scheme::OneStepScheme, rigid::Symbol, callback, coarsen::Bool,
                                         verbose::Bool, ctrl)
    @nospecialize curve U dt symbol saveat scheme callback ctrl
    several = curve isa AbstractVector
    problem = CurveProblem(several ? curve : [curve], _SimpleVelocity(U, several); symbol=symbol, rigid=rigid)
    integ = init(problem, scheme; t0=t0, tol=tol, Nmax=Nmax)
    _drive!(integ, t1, dt, saveat, tol, Nmax, callback, coarsen, several, verbose, ctrl)
end

"""
    solve(problem, scheme = ETDRK2(); tspan, dt, saveat = nothing, tol = 1e-6, Nmax = 16384, callback = nothing,
          coarsen = false, verbose = false, adaptive = false, errtol = 1e-5, dtmin = 0.0, dtmax = Inf) -> Evolution

Run `problem` over `tspan` as `evolve` does (`dt` a number or `dt(Xs, t)`; or `adaptive = true`, with `dt`
optional), from `init(problem, scheme)`.
"""
function solve(problem::CurveProblem, scheme::OneStepScheme=ETDRK2(); tspan, dt=nothing, saveat=nothing,
               tol::Real=1e-6, Nmax::Int=16384, callback=nothing, coarsen::Bool=false, verbose::Bool=false,
               adaptive::Bool=false, errtol::Real=1e-5, dtmin::Real=0.0, dtmax::Real=Inf)
    ctrl = _step_control(dt, scheme, adaptive, errtol, dtmin, dtmax)
    integ = init(problem, scheme; t0=tspan[1], tol=tol, Nmax=Nmax)
    _drive!(integ, Float64(tspan[2]), dt, saveat, Float64(tol), Nmax, callback, coarsen, true, verbose, ctrl)
end

# The step controller of `adaptive` (`nothing` without it): its bounds, the step it would take next (`h`, NaN until
# chosen) and whether the last attempt was rejected (no growth right after a rejection)
mutable struct _StepControl
    errtol::Float64
    dtmin::Float64
    dtmax::Float64
    h::Float64
    rejected::Bool
end
const _SAFETY = 0.9          # the next step aims at 0.9² ≈ 0.8 errtol: few rejections, little lost
const _FMIN = 0.2            # the step shrinks at most 5× after a step,
const _FMAX = 5.0            # and grows at most 5× (out of an initial layer in a few steps)
const _FFAIL = 0.25          # after a rejection that gives no estimate (StepFailure, a non-positive length)
const _MAX_REJECTIONS = 30   # in a row: 0.25³⁰ ≈ 1e-18, the step is then not the trouble
function _step_control(dt, scheme::OneStepScheme, adaptive::Bool, errtol::Real, dtmin::Real, dtmax::Real)
    if !adaptive
        dt === nothing && throw(ArgumentError("give the step dt (a number or dt(X, t)), or adaptive = true"))
        return nothing
    end
    has_error_estimate(scheme) || throw(ArgumentError("adaptive = true needs a scheme with an error estimate " *
                                                      "(ETDRK2), not " * string(nameof(typeof(scheme)))))
    errtol > 0 || throw(ArgumentError("errtol = " * string(errtol) * ": must be positive"))
    0 <= dtmin < dtmax || throw(ArgumentError("dtmin = " * string(dtmin) * ", dtmax = " * string(dtmax) *
                                              ": need 0 ≤ dtmin < dtmax"))
    dt isa Number && !(dtmin <= dt <= dtmax) &&
        throw(ArgumentError("the first step dt = " * string(dt) * " is not in [dtmin, dtmax]"))
    _StepControl(Float64(errtol), Float64(dtmin), Float64(dtmax), dt isa Number ? Float64(dt) : NaN, false)
end

# The first step when none is given: the time in which the fastest point moves sqrt(errtol)/10 of its curve's
# length (the estimate is ~ (h·rate)² for a smooth motion); a step too large is rejected at the cost of one
# evaluation and shrunk
function _first_step(integ::ThetaLIntegrator, errtol::Float64)
    integ.U === nothing && _evaluate_initial!(integ)
    rate = 0.0
    for (X, U) in zip(integ.states, integ.U::Vector{Vector{Float64}})
        isempty(U) || (rate = max(rate, maximum(abs, U) / X.L))
    end
    rate > 0 ? 0.1 * sqrt(errtol) / rate : Inf
end

# The loop of `evolve` and `solve` (functions positional and unspecialized: one compiled loop for all of them)
Base.@nospecializeinfer function _drive!(integ::ThetaLIntegrator, t1::Float64, dt, saveat, tol::Float64, Nmax::Int,
                                         callback, coarsen::Bool, several::Bool, verbose::Bool, ctrl=nothing)
    @nospecialize dt saveat callback ctrl
    keep() = several ? copy(integ.states) : integ.states[1]
    wall = time(); last = (integ.t, integ.stats.steps)        # the time and steps at the previous line
    function report()
        verbose || return
        n = integ.stats.steps - last[2]; dtmean = n > 0 ? (integ.t - last[1]) / n : NaN
        println(progress_line(integ; dt=dtmean) * "  " * string(round(time() - wall; digits=1)) * " s")
        flush(stdout)                                  # seen at once when the output goes to a file
        last = (integ.t, integ.stats.steps)
    end
    report()
    ts = [integ.t]; out = Any[keep()]; ids = [copy(integ.ids)]
    fs = NamedTuple[integ.fields]; ys = Union{Nothing,Vector{Float64}}[integ.bulk]   # (replaced, never changed)
    calm = Dict{Int,Int}()                             # per id: steps the curve could have been halved
    stopped = false
    for stop in _stops(saveat, integ.t, t1)
        while integ.t < stop
            t = integ.t
            if ctrl === nothing
                h = min(dt isa Number ? Float64(dt) : Float64(dt(several ? integ.states : integ.states[1], t)),
                        stop - t)
                stop - t - h <= 1e-9 * h && (h = stop - t)                # no sliver of a step left by rounding
                status = advance!(integ, h == stop - t ? stop : t + h)
                status === :success || error("the step from t = $t was rejected ($(integ.stats.last_rejection)); " *
                                             "evolve and solve do not retry with a given dt: use adaptive = " *
                                             "true, or init and advance! to control steps")
            else
                _adaptive_step!(integ, ctrl::_StepControl, stop, dt, several) || continue   # rejected: retry
            end
            _adapt_resolution!(integ, tol, Nmax, coarsen, calm, verbose)
            if callback !== nothing && callback(integ) === :stop
                stopped = true; break
            end
        end
        push!(ts, integ.t); push!(out, keep()); push!(ids, copy(integ.ids)); push!(fs, integ.fields)
        push!(ys, integ.bulk)
        report()
        stopped && break
    end
    st = integ.stats
    several ? Evolution(ts, Vector{Vector{ThetaLState}}(out), ids, st.steps, st.rejections, fs, ys) :
              Evolution(ts, Vector{ThetaLState}(out), ids, st.steps, st.rejections, fs, ys)
end

# One attempt of the adaptive loop towards `stop`: true if accepted. The step taken is the controller's, bounded by
# a step function dt, and fitted to the time left: all of it if within 10% of the step (and the bound), half of it
# if less than two steps (no sliver before a saved time). After the attempt the controller's step is the next one: from the estimate
# (err ∝ h² for the embedded order-1 step, hence the square root; err = 0 grows it the most), back up to the step it
# had if the attempt was cut short by a saved time and the estimate allows, never grown right after a rejection.
Base.@nospecializeinfer function _adaptive_step!(integ::ThetaLIntegrator, ctrl::_StepControl, stop::Float64, dt,
                                                 several::Bool)
    @nospecialize dt
    t = integ.t; left = stop - t
    isnan(ctrl.h) && (ctrl.h = clamp(_first_step(integ, ctrl.errtol), ctrl.dtmin, ctrl.dtmax))
    hmax = ctrl.dtmax
    dt isa Function && (hmax = max(ctrl.dtmin, min(hmax, Float64(dt(several ? integ.states : integ.states[1], t)))))
    h = min(ctrl.h, hmax)
    hc = left <= min(1.1h, hmax) ? left : left < 2h ? left / 2 : h          # (never stretched past the bound)
    status = advance!(integ, hc == left ? stop : t + hc; errtol=ctrl.errtol)
    st = integ.stats; err = st.error_estimate
    if status === :success
        ρ = err > 0 ? _SAFETY * sqrt(ctrl.errtol / err) : isnan(err) ? 1.0 : Inf
        ctrl.h = clamp(hc * clamp(ρ, _FMIN, ctrl.rejected ? 1.0 : max(_FMAX, h / hc)), ctrl.dtmin, ctrl.dtmax)
        ctrl.rejected = false
        return true
    end
    by_estimate = st.last_rejection === :error_estimate && err > 0 && isfinite(err)
    ctrl.h = hc * (by_estimate ? clamp(_SAFETY * sqrt(ctrl.errtol / err), _FMIN, _SAFETY) : _FFAIL)
    ctrl.rejected = true
    why = "the step from t = " * string(t) * " was rejected " * string(st.consecutive_rejections) *
          " times in a row (last: " * string(st.last_rejection) *
          (by_estimate ? ", error estimate " * _short(err, 3) * " > errtol" : "") * ")"
    st.consecutive_rejections >= _MAX_REJECTIONS && error(why)
    (ctrl.h < ctrl.dtmin || t + ctrl.h == t) &&
        error(why * ", and the next step " * _short(ctrl.h, 3) * " would be below dtmin = " * string(ctrl.dtmin))
    false
end

# N doubled where the tail exceeds tol; with `coarsen`, halved where the tail on half the points has stayed below
# tol/100 for 10 steps (the gap between the two thresholds and the wait keep them from undoing each other);
# with `verbose`, each change is printed
function _adapt_resolution!(integ::ThetaLIntegrator, tol::Float64, Nmax::Int, coarsen::Bool, calm::Dict{Int,Int},
                            verbose::Bool)
    for c in eachindex(integ.states)
        X = integ.states[c]; N = length(X); id = integ.ids[c]; tail = spectral_tail(X)
        if tail > tol && 2N <= Nmax
            double_grid!(integ, c); calm[id] = 0
            verbose && _report_resolution(integ, id, N, 2N, "tail " * _short(tail, 2) * " > tol")
        elseif coarsen && N >= 128 && spectral_tail(halve_grid(X)) < tol / 100
            calm[id] = get(calm, id, 0) + 1
            if calm[id] >= 10
                halve_grid!(integ, c); calm[id] = 0
                verbose && _report_resolution(integ, id, N, N ÷ 2, "smooth for 10 steps")
            end
        else
            calm[id] = 0
        end
    end
    integ
end
function _report_resolution(integ::ThetaLIntegrator, id, N::Int, Nnew::Int, why::String)
    curve = length(integ.states) == 1 ? "" : "curve " * string(id) * ": "
    println("t = ", _short(integ.t, 4), "  ", curve, "N ", N, " → ", Nnew, " (", why, ")")
    flush(stdout)
end

"""
    save_curves(path_or_io, sol; times=sol.t)

Write the curves of `sol` at the saved `times` as CSV: a header `t,curve,x,y`, then one row per point (`curve` is
the curve's id), counterclockwise from x₀; the last point is not repeated.
"""
function save_curves(io::IO, sol::Evolution; times=sol.t)
    println(io, "t,curve,x,y")
    for t in times
        i = findfirst(τ -> abs(τ - t) <= 1e-12 * (1 + max(abs(τ), abs(t))), sol.t)
        i === nothing && throw(ArgumentError("no curves were saved at t = $t"))
        Xs = sol.curves[i]
        for (j, X) in enumerate(Xs isa ThetaLState ? [Xs] : Xs), x in curve_points(X)
            println(io, t, ",", sol.ids[i][j], ",", x[1], ",", x[2])
        end
    end
end
save_curves(path::AbstractString, sol::Evolution; times=sol.t) = (open(io -> save_curves(io, sol; times), path, "w"); path)
