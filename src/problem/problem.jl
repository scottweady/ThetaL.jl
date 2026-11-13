# The problem: N closed curves moving with normal velocities that depend on all of them,
# (U₁, …, U_N) = F(X₁, …, X_N; t), with a stiff symbol per curve, and the state carried with them (fields on the
# curves, a bulk state; problem/carried.jl).

"""
    CurveProblem(curves, velocity; symbol = nothing, rigid = :full, tags = nothing, validate = nothing,
                 fields = (;), bulk = nothing)

What is evolved:

- `curves`: a vector of closed counterclockwise curves, each a parametrization `ξ -> (x, y)` (2π-periodic) or a
  `ThetaLState`;
- `velocity(states, ev::Evaluation; warm) -> U` or `(; U, ctx, W, rates, ydot)`: the outward normal velocity of
  every curve at the grid points of its state (one vector per curve), all curves at once; `ctx` is any context the
  caller wants back as `warm` at the next evaluation (a warm start). With carried state: `W` the material tangential
  velocity (one vector per curve; needed when a field is transported), `rates` the fields' rates (a named tuple,
  one vector per curve each; a field left out has rate 0) and `ydot` the bulk state's rate. Throw
  `StepFailure(reason)` to reject a step;
- `symbol`: the growth rate σ of a small-scale mode of a curve, θ_t ≈ σ φ̂, as a function `σ(q)` of the wavenumber
  q = 2π|k|/L along the curve, or `σ(q, X)` of it and the curve X (a `ThetaLState`), frozen per curve at the start
  of each step: `q -> -q^2` for U = −κ; `nothing` for none;
- `rigid`: which rigid motion is taken out of U and integrated exactly, `:full` (translation and rotation),
  `:translation` or `:none`;
- `tags`: one opaque tag per curve, carried through edits and checkpoints, never interpreted;
- `validate(states) -> Bool`: optional; a step whose end state fails it is rejected;
- `fields`: fields on the curves, a named tuple of `CurveField`s (`(; c = CurveField(c0))`); their values at each
  stage are `ev.fields.c` (one vector per curve), at the integrator's time `integ.fields.c`;
- `bulk`: a `BulkState`, one vector carried with the curves; its stage value is `ev.bulk`, `integ.bulk` its value.

`velocity` and `symbol` may be replaced between steps (`prob.velocity = …`); the change applies from the next step.
"""
mutable struct CurveProblem
    curves::Vector{Any}
    velocity::Any
    symbol::Any
    rigid::Symbol
    tags::Vector{Any}
    validate::Any
    field_names::Vector{Symbol}
    fields::Vector{CurveField}
    bulk::Union{Nothing,BulkState}
end
Base.@nospecializeinfer function CurveProblem(curves::AbstractVector, velocity; symbol=nothing, rigid::Symbol=:full,
                                              tags=nothing, validate=nothing, fields::NamedTuple=(;),
                                              bulk::Union{Nothing,BulkState}=nothing)
    @nospecialize curves velocity symbol tags validate fields
    rigid in (:full, :translation, :none) || throw(ArgumentError("rigid must be :full, :translation or :none"))
    t = tags === nothing ? Any[nothing for _ in curves] : Vector{Any}(collect(tags))
    length(t) == length(curves) || throw(DimensionMismatch("$(length(t)) tags for $(length(curves)) curves"))
    names = Symbol[k for k in keys(fields)]; fs = CurveField[]
    for k in names
        f = getproperty(fields, k)
        f isa CurveField || throw(ArgumentError("field $k must be a CurveField"))
        push!(fs, f)
    end
    CurveProblem(Vector{Any}(collect(curves)), velocity, symbol, rigid, t, validate, names, fs, bulk)
end
