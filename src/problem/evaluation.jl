# The velocity protocol: what the caller's velocity function is told (`Evaluation`), how it signals a failure
# (`StepFailure`), and how its result is read (`_evaluate_velocity`).
#
#   velocity(states, ev::Evaluation; warm) -> U                    (one Vector per curve)
#   velocity(states, ev::Evaluation; warm) -> (; U, ctx)           (ctx is handed back as `warm` next time)
#   velocity(states, ev::Evaluation; warm) -> (; U, ctx, W, rates, ydot)    (with state carried with the curves)
#
# `Evaluation` may gain fields; a named-tuple result may gain fields. Callers that do not use them are not affected.

"""
    Evaluation

What a velocity evaluation is for, passed to `velocity(states, ev; warm)`:

- `t`: the time of the evaluation (exact: at the end of a step, the time stepped to);
- `kind`: `:init` (the start of a run, or the first evaluation after a curve edit), `:stage` (inside a step) or
  `:end` (the end of a step; its result is also the next step's first stage); `:probe` for a perturbed state of
  `measure_symbol`/`check_symbol` (no step: `stage = 0`, `h = 0.0`; nothing of it is kept);
- `stage`, `c`: the stage index and its fraction of the step (`0` and `0.0` at `:init`, `c = 1` at `:end`);
- `h`: the step being taken (`0.0` at `:init`);
- `attempt`: 1, or more when the step is a retry after a rejection;
- `ids`, `tags`: per curve, a stable id and the caller's opaque tag (ThetaL never interprets tags);
- `topology`: a counter bumped by every curve edit and change of N; a warm start keyed on it knows when the
  curves have changed;
- `base`: the states at the start of the step (a moving bulk mesh forms its displacement from them);
- `fields`: the fields on the curves at this stage, a named tuple (`ev.fields.c`: one vector per curve);
- `bulk`: the bulk state at this stage (`nothing` without one).
"""
struct Evaluation
    t::Float64
    kind::Symbol
    stage::Int
    c::Float64
    h::Float64
    attempt::Int
    ids::Vector{Int}
    tags::Vector{Any}
    topology::Int
    base::Vector{ThetaLState}
    fields::NamedTuple
    bulk::Union{Nothing,Vector{Float64}}
end

"""
    StepFailure(reason)

Thrown by a velocity function to reject the step it is evaluated for (a solver that did not converge, a state it
cannot handle). `advance!` catches it, leaves the integrator as it was and returns `:rejected`, with `reason` in
`integ.stats.last_rejection`. Any other exception propagates (the integrator is still left as it was).
"""
struct StepFailure <: Exception
    reason::Any
end
Base.showerror(io::IO, e::StepFailure) = print(io, "StepFailure: ", e.reason)

# Call the problem's velocity and read its result: (U::Vector{Vector{Float64}}, ctx, aux::_Aux). U is copied (so a
# caller that reuses its buffers cannot change the integrator's data) and checked against the curves' sizes; aux
# holds W, the fields' rates and ydot when the problem carries state (problem/carried.jl).
Base.@nospecializeinfer function _evaluate_velocity(problem, states::Vector{ThetaLState}, ev::Evaluation, warm)
    @nospecialize warm
    r = problem.velocity(states, ev; warm=warm)
    if r isa NamedTuple
        haskey(r, :U) || throw(ArgumentError("the velocity returned a named tuple without U"))
        Ur = r.U; ctx = haskey(r, :ctx) ? r.ctx : nothing
    elseif r isa Tuple
        throw(ArgumentError("the velocity returned a Tuple: return U (one vector per curve) or a named tuple " *
                            "(; U, ctx); the tangential velocity is not used"))
    else
        Ur = r; ctx = nothing
    end
    carried = !isempty(problem.fields) || problem.bulk !== nothing      # (behind a barrier: compiled only if used)
    _checked_velocities(Ur, states), ctx,
        carried ? Base.inferencebarrier(_read_aux)(problem, r, states)::_Aux : _NO_AUX
end

Base.@nospecializeinfer function _checked_velocities(Ur, states::Vector{ThetaLState})
    @nospecialize Ur
    length(Ur) == length(states) ||
        throw(DimensionMismatch("the velocity returned $(length(Ur)) vectors for $(length(states)) curves"))
    U = Vector{Vector{Float64}}(undef, length(states))
    for c in eachindex(states)
        u = Vector{Float64}(Ur[c])                         # a copy
        length(u) == length(states[c].φ) ||
            throw(DimensionMismatch("curve $c has $(length(states[c].φ)) points, its velocity $(length(u)) values"))
        U[c] = u
    end
    U
end
