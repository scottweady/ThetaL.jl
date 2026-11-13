# The drivers: run an integrator over a time span with a fixed step (or a step function), landing on the saved
# times, adapting N, and keeping the curves at the saved times. `evolve` for a curve (or several) and a simple
# velocity U(X, t); `solve` for a `CurveProblem`. Their result (`Evolution`, `sol(t)`) and a CSV writer for it
# (`save_curves`).

"""
    Evolution

The result of `evolve` and `solve`: the saved times `t`, the curves there (`curves`: a `ThetaLState` each for
`evolve` of one curve, otherwise a vector of states per time), the curve `ids` at each saved time, and the number of
accepted `steps`; with state carried with the curves, `fields` (per saved time a named tuple, one vector per curve
each: `sol.fields[i].c`) and `bulk` (per saved time). `sol(t)` is the curves saved at time `t`. With DotPlots
loaded, `dotplot(sol)` draws them in the terminal, one frame per saved time.
"""
struct Evolution{C}
    t::Vector{Float64}
    curves::Vector{C}
    ids::Vector{Vector{Int}}
    steps::Int
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
           rigid = :full, callback = nothing, coarsen = false, verbose = false) -> Evolution

Move the closed counterclockwise `curve` over `tspan = (t0, t1)` with outward normal velocity `U`. One call for
simple scripts; `init`/`advance!` for anything more.

- `curve`: a parametrization `ξ -> (x, y)` (2π-periodic), a `ThetaLState`, or a vector of either: several curves
  that move together.
- `U(X, t)`: the normal velocity at the N points of the curve `X` (a `ThetaLState`: use `curve_points(X)`,
  `curvature(X)`, `tangent_angles(X)`, …; for a velocity field u, `velocity_components(X, u)[1]`); for several
  curves `U(Xs, t)`, one vector per curve. Its rigid part is integrated exactly (`rigid`, as in `CurveProblem`).
- `symbol`: the growth rate σ of a small-scale mode, θ_t ≈ σ φ̂, as `σ(q)` of the wavenumber q = 2π|k|/L along
  the curve or `σ(q, X)`: `q -> -q^2` for U = −κ, `q -> -q^4` for U = κ_ss; `nothing`: none.
- `dt`: the step, a number or a function `dt(X, t)` (`dt(Xs, t)` for several curves).
- `saveat`: the times to keep: `nothing` (the start and the end), a spacing, or a list. Steps land on them exactly.
- `tol`: N is the smallest power of two from 64 whose φ tail (`spectral_tail`) is below tol/100; it is doubled (up
  to `Nmax`) after any step whose tail exceeds `tol`, and with `coarsen = true` halved when the tail on half the
  points has stayed below tol/100 for 10 steps. Without `coarsen`, N never comes back down: one step of noise in a
  computed velocity is enough to double it for the rest of the run (`ETDRK2(filter_cutoff = :auto)` damps that
  noise). Set `Nmax` to what your velocity can afford: for an O(N²) solve the default 16384 is out of reach.
- `callback(integ)`: called after every step; return `:stop` to end the run. It may edit the curves.
- `verbose`: print a line at the start and at each saved time (`progress_line`: t, steps, the mean Δt since the
  previous line, N, evaluations), with the wall time so far, and a line whenever N changes (when, and why).

A rejected step is an error here (`evolve` does not retry); use the integrator to control steps.
"""
function evolve(curve, U, tspan; dt, symbol=nothing, saveat=nothing, tol::Real=1e-6, Nmax::Int=16384,
                scheme::OneStepScheme=ETDRK2(), rigid::Symbol=:full, callback=nothing, coarsen::Bool=false,
                verbose::Bool=false)
    _evolve(curve, U, Float64(tspan[1]), Float64(tspan[2]), dt, symbol, saveat, Float64(tol), Nmax, scheme, rigid,
            callback, coarsen, verbose)
end
# (the caller's functions positional and unspecialized: one compiled `evolve` for every script, whose own closures
# cost only the keyword wrapper above)
Base.@nospecializeinfer function _evolve(curve, U, t0::Float64, t1::Float64, dt, symbol, saveat, tol::Float64,
                                         Nmax::Int, scheme::OneStepScheme, rigid::Symbol, callback, coarsen::Bool,
                                         verbose::Bool)
    @nospecialize curve U dt symbol saveat scheme callback
    several = curve isa AbstractVector
    problem = CurveProblem(several ? curve : [curve], _SimpleVelocity(U, several); symbol=symbol, rigid=rigid)
    integ = init(problem, scheme; t0=t0, tol=tol, Nmax=Nmax)
    _drive!(integ, t1, dt, saveat, tol, Nmax, callback, coarsen, several, verbose)
end

"""
    solve(problem, scheme = ETDRK2(); tspan, dt, saveat = nothing, tol = 1e-6, Nmax = 16384, callback = nothing,
          coarsen = false, verbose = false) -> Evolution

Run `problem` over `tspan` as `evolve` does (`dt` a number or `dt(Xs, t)`), from `init(problem, scheme)`.
"""
function solve(problem::CurveProblem, scheme::OneStepScheme=ETDRK2(); tspan, dt, saveat=nothing,
               tol::Real=1e-6, Nmax::Int=16384, callback=nothing, coarsen::Bool=false, verbose::Bool=false)
    integ = init(problem, scheme; t0=tspan[1], tol=tol, Nmax=Nmax)
    _drive!(integ, Float64(tspan[2]), dt, saveat, Float64(tol), Nmax, callback, coarsen, true, verbose)
end

# The loop of `evolve` and `solve` (functions positional and unspecialized: one compiled loop for all of them)
Base.@nospecializeinfer function _drive!(integ::ThetaLIntegrator, t1::Float64, dt, saveat, tol::Float64, Nmax::Int,
                                         callback, coarsen::Bool, several::Bool, verbose::Bool)
    @nospecialize dt saveat callback
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
            h = min(dt isa Number ? Float64(dt) : Float64(dt(several ? integ.states : integ.states[1], t)), stop - t)
            stop - t - h <= 1e-9 * h && (h = stop - t)                # no sliver of a step left by rounding
            status = advance!(integ, h == stop - t ? stop : t + h)
            status === :success || error("the step from t = $t was rejected ($(integ.stats.last_rejection)); " *
                                         "evolve and solve do not retry: use init and advance! to control steps")
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
    several ? Evolution(ts, Vector{Vector{ThetaLState}}(out), ids, integ.stats.steps, fs, ys) :
              Evolution(ts, Vector{ThetaLState}(out), ids, integ.stats.steps, fs, ys)
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
