# The integrator: the state of a run (time, curves, their velocities and the caller's context, ids and tags),
# advanced one step at a time with exact landing (`advance!`), transactional (a rejected step changes nothing).
# Curve edits are in edits.jl, checkpoints in checkpoint.jl, the drivers (`evolve`, `solve`) in src/drivers/.

"""Counts of a run: `steps` (accepted), `evaluations` (velocity calls), `rejections`, `consecutive_rejections`,
`error_estimate` of the last step tried (`NaN` if the scheme has none), `last_rejection` (its reason), `last_step`
(the size of the last accepted step; `NaN` before the first)."""
mutable struct IntegratorStats
    steps::Int
    evaluations::Int
    rejections::Int
    consecutive_rejections::Int
    error_estimate::Float64
    last_rejection::Any
    last_step::Float64
end
IntegratorStats() = IntegratorStats(0, 0, 0, 0, NaN, nothing, NaN)

"""
    ThetaLIntegrator

The state of a run, made by `init(problem, scheme)`. Read-only fields (change them only through `advance!`,
`step!` and the curve edits):

- `t`: the time; `states::Vector{ThetaLState}`: the curves; `U`: their normal velocities at `t` (`nothing` after
  an edit, until the next step evaluates them);
- `ctx`: the context of the last committed velocity evaluation (handed to the next as `warm`);
- `ids`, `tags`: per curve, a stable id and the caller's tag; `topology`: bumped by every edit and change of N;
- `rigid::Vector{RigidMotion}`: each curve's rigid motion over the last step;
- `fields`: the fields on the curves at `t`, a named tuple (`integ.fields.c`: one vector per curve); `bulk`: the
  bulk state at `t` (`nothing` without one);
- `stats::IntegratorStats`.

Settable between steps: `scheme` (or its filter), `problem.velocity`, `problem.symbol`; they apply from the next
step. Vectors are replaced, never changed in place, so a reference to `integ.states` or `integ.U` stays valid.
"""
mutable struct ThetaLIntegrator
    problem::CurveProblem
    scheme::OneStepScheme
    t::Float64
    states::Vector{ThetaLState}
    U::Union{Nothing,Vector{Vector{Float64}}}
    ctx::Any
    ids::Vector{Int}
    tags::Vector{Any}
    topology::Int
    next_id::Int
    rigid::Vector{RigidMotion}
    stats::IntegratorStats
    tol::Float64                                    # resolution of curves made from parametrizations
    Nmax::Int
    fields::NamedTuple                              # the fields on the curves (one vector per curve each)
    bulk::Union{Nothing,Vector{Float64}}
    aux::_Aux                                       # W, rates and ydot of the evaluation that gave U
end

# the fields as vectors, fields[f][c], in the problem's order (behind a barrier, as in the step: compiled only for
# problems with fields)
_field_vectors(integ::ThetaLIntegrator) = isempty(integ.problem.field_names) ? Vector{Vector{Float64}}[] :
    Base.inferencebarrier(_field_vectors_)(integ.fields, integ.problem.field_names)::Vector{Vector{Vector{Float64}}}
_field_vectors_(fields::NamedTuple, names::Vector{Symbol}) =
    Vector{Vector{Float64}}[getproperty(fields, k) for k in names]
function _set_fields!(integ::ThetaLIntegrator, fv::Vector{Vector{Vector{Float64}}})
    isempty(fv) && return integ
    integ.fields = Base.inferencebarrier(_fields_tuple)(integ.problem.field_names, fv)::NamedTuple
    integ
end

"""
    init(problem, scheme = ETDRK2(); t0 = 0.0, U = nothing, ctx = nothing, from = nothing, tol = 1e-6,
         Nmax = 16384) -> ThetaLIntegrator

Start a run of `problem` with `scheme`. Curves given as parametrizations are resolved to `tol` (N a power of two
from 64, at most `Nmax`). The velocities at `t0` are evaluated (`kind = :init`, `warm = ctx`) unless `U` is given.
`from = checkpoint(integ)` restarts a run instead (its time, curves, velocities, ids, tags, resolved filter cutoff
and carried state). The fields start from the problem's `CurveField`s, the bulk state from its `BulkState`.
"""
function init(problem::CurveProblem, scheme::OneStepScheme=ETDRK2(); t0::Real=0.0, U=nothing,
                          ctx=nothing, from=nothing, tol::Real=1e-6, Nmax::Int=16384)
    from === nothing || return _init_from_checkpoint(problem, scheme, from, Float64(tol), Nmax)
    states = ThetaLState[c isa ThetaLState ? c : _from_curve(c, Float64(tol), 64, Nmax) for c in problem.curves]
    n = length(states)
    integ = ThetaLIntegrator(problem, scheme, Float64(t0), states, nothing, ctx, collect(1:n), copy(problem.tags),
                             0, n + 1, fill(NO_RIGID_MOTION, n), IntegratorStats(), Float64(tol), Nmax,
                             isempty(problem.fields) ? _NO_FIELDS :
                             Base.inferencebarrier(_initial_fields_tuple)(problem, states)::NamedTuple,
                             problem.bulk === nothing ? nothing : copy(problem.bulk.y0), _NO_AUX)
    U === nothing || (problem.bulk === nothing && isempty(problem.fields)) ||
        throw(ArgumentError("with state carried with the curves, give no U: the first evaluation also gives W, " *
                            "the rates and ydot"))
    U === nothing ? _evaluate_initial!(integ) : (integ.U = _checked_velocities(U, states))
    _resolve_auto_cutoff!(integ)
    integ
end

# the velocities at the integrator's time (kind :init): at the start, and after a curve edit
function _evaluate_initial!(integ::ThetaLIntegrator)
    ev = Evaluation(integ.t, :init, 0, 0.0, 0.0, 1, integ.ids, integ.tags, integ.topology, integ.states,
                    integ.fields, integ.bulk)
    integ.U, integ.ctx, integ.aux = _evaluate_velocity(integ.problem, integ.states, ev, integ.ctx)
    integ.stats.evaluations += 1
    integ
end

# a filter cutoff `:auto` becomes an integer from the noise in the first velocities (`noise_cutoff`)
function _resolve_auto_cutoff!(integ::ThetaLIntegrator)
    f = integ.scheme.filter
    if f.cutoff === :auto && integ.U !== nothing
        kc = noise_cutoff(integ.states, integ.U, integ.problem.symbol)
        integ.scheme = with_filter(integ.scheme, SpectralFilter(order=f.order, cutoff=kc))
    end
    integ
end

"""
    advance!(integ, t; errtol = Inf) -> :success or :rejected

One step from `integ.t` to `t` (exactly: `integ.t == t` afterwards). If the velocities are not current (after a
curve edit) they are evaluated first (`kind = :init`). A step is rejected, and the integrator left exactly as it
was, when the velocity throws `StepFailure`, a stage has a non-positive length, the embedded error estimate exceeds
`errtol` (checked before the end evaluation, so that costs one evaluation less), or `validate` fails; the reason is
in `integ.stats.last_rejection`, and a retry is told `attempt = 2, 3, …`.
"""
function advance!(integ::ThetaLIntegrator, t::Real; errtol::Real=Inf)
    t_end = Float64(t)
    t_end > integ.t || throw(ArgumentError("advance! to t = $t_end: not after the integrator's time $(integ.t)"))
    integ.U === nothing && _evaluate_initial!(integ)
    _resolve_auto_cutoff!(integ)
    st = integ.stats
    r = _thetaL_step(integ.scheme, integ.problem, integ.states, integ.U::Vector{Vector{Float64}}, integ.t, t_end,
                     integ.ctx, integ.ids, integ.tags, integ.topology, st.consecutive_rejections + 1, Float64(errtol),
                     _field_vectors(integ), integ.bulk, integ.aux)
    st.evaluations += r.evaluations; st.error_estimate = r.error_estimate
    if r.status === :success
        st.last_step = t_end - integ.t
        integ.t = t_end; integ.states = r.states; integ.U = r.U; integ.ctx = r.ctx; integ.rigid = r.rigid
        _set_fields!(integ, r.fields); integ.bulk = r.bulk; integ.aux = r.aux
        st.steps += 1; st.consecutive_rejections = 0
    else
        st.rejections += 1; st.consecutive_rejections += 1; st.last_rejection = r.reason
    end
    r.status
end

"""
    step!(integ, dt; errtol = Inf) -> :success or :rejected

`advance!(integ, integ.t + dt; errtol)`. Use `advance!` to land exactly on a given time.
"""
step!(integ::ThetaLIntegrator, dt::Real; errtol::Real=Inf) = advance!(integ, integ.t + dt; errtol=errtol)

Base.show(io::IO, integ::ThetaLIntegrator) =
    print(io, "ThetaLIntegrator(", nameof(typeof(integ.scheme)), ", t = ", integ.t, ", ", length(integ.states),
          " curve", length(integ.states) == 1 ? "" : "s", ", steps = ", integ.stats.steps, ")")

# x to `digits` significant digits, as text (no Printf: its formats cost a second of compilation)
_short(x::Real, digits::Int) = string(round(Float64(x); sigdigits=digits))

"""
    progress_line(integ; dt = integ.stats.last_step) -> String

One line on the run so far: the time, the accepted steps (and the rejected ones, if any), the step Δt (by default
the last one), the number of points N (its range and the number of curves, for several) and the velocity
evaluations. Printed by `evolve` and `solve` with `verbose = true` (with Δt the mean step since the previous line);
print it yourself when driving the integrator by hand.
"""
function progress_line(integ::ThetaLIntegrator; dt::Real=integ.stats.last_step)
    st = integ.stats; n = length(integ.states); Nmin = typemax(Int); Nmax = 0
    for X in integ.states
        Nmin = min(Nmin, length(X.φ)); Nmax = max(Nmax, length(X.φ))
    end
    # (strings joined, not one `string` of mixed arguments: that compiles a method per combination of types)
    N = n == 0 ? "no curves" : n == 1 ? "N " * string(Nmax) :
        "N " * string(Nmin) * "–" * string(Nmax) * " (" * string(n) * " curves)"
    rej = st.rejections > 0 ? " (" * string(st.rejections) * " rejected)" : ""
    Δt = isnan(dt) ? "—" : _short(dt, 3)
    join(String["t = ", rpad(_short(integ.t, 4), 9), "  steps ", rpad(string(st.steps), 6), rej, "  Δt ", rpad(Δt, 8),
                "  ", rpad(N, 8), "  evaluations ", string(st.evaluations)])
end
