# Edits of the integrator's curves between steps. Because every curve's velocity may depend on every curve
# (a coupled F), any change of a curve's geometry or of the set of curves invalidates the stored velocities of all
# curves: the next `advance!` evaluates them first (`kind = :init`). Changing N (`double_grid!`, `halve_grid!`)
# keeps the curve, so its velocity is interpolated instead. Every edit bumps `topology`. Fields on the curves
# follow their curve: interpolated with it, resampled when it is replaced by one of another N, given (or taken from
# the field's `init`) for an added curve, dropped with a removed one.
#
# Curves are addressed by their position `c` in `integ.states`; `curve_index(integ, id)` finds a curve's position
# from its stable id.

"""The position in `integ.states` of the curve with id `id` (`nothing` if there is none)."""
curve_index(integ::ThetaLIntegrator, id::Int) = findfirst(==(id), integ.ids)

# replace curve c (vectors replaced, not changed in place)
function _replace_curve(integ::ThetaLIntegrator, c::Int, X::ThetaLState)
    states = copy(integ.states); states[c] = X; integ.states = states
end
function _invalidate_velocities!(integ::ThetaLIntegrator)
    integ.U = nothing; integ.aux = _NO_AUX; integ.topology += 1
    integ
end

# g applied to curve c's entry of each field, and of W and the rates (when current)
function _map_curve!(g, integ::ThetaLIntegrator, c::Int)
    fv = _field_vectors(integ)
    for f in eachindex(fv)
        v = copy(fv[f]); v[c] = g(v[c]); fv[f] = v
    end
    _set_fields!(integ, fv)
    a = integ.aux
    if integ.U !== nothing && a !== _NO_AUX
        W = a.W === nothing ? nothing : (v = copy(a.W); v[c] = g(v[c]); v)
        rates = Union{Nothing,Vector{Vector{Float64}}}[r === nothing ? nothing : (v = copy(r); v[c] = g(v[c]); v)
                                                       for r in a.rates]
        integ.aux = _Aux(W, rates, a.ydot)
    end
    integ
end

"""Curve `c` on twice the points (its velocity interpolated with it)."""
function double_grid!(integ::ThetaLIntegrator, c::Int)
    _replace_curve(integ, c, double_grid(integ.states[c]))
    integ.U === nothing || (U = copy(integ.U); U[c] = double_grid(U[c]); integ.U = U)
    _map_curve!(double_grid, integ, c)
    integ.topology += 1
    integ
end

"""Curve `c` on half the points (spectral truncation; its velocity truncated with it)."""
function halve_grid!(integ::ThetaLIntegrator, c::Int)
    _replace_curve(integ, c, halve_grid(integ.states[c]))
    integ.U === nothing || (U = copy(integ.U); U[c] = halve_grid(U[c]); integ.U = U)
    _map_curve!(halve_grid, integ, c)
    integ.topology += 1
    integ
end

"""
    set_curve!(integ, c, curve; fields = nothing)

Replace curve `c` by `curve` (a `ThetaLState` or a parametrization); its id and tag are kept. Its fields are kept,
resampled to the new N if it changes, or set from `fields` (a named tuple: per field a number, a vector or a
function `X -> values`).
"""
function set_curve!(integ::ThetaLIntegrator, c::Int, curve; fields=nothing)
    X = _as_state(integ, curve)
    _replace_curve(integ, c, X)
    fv = _field_vectors(integ)
    for f in eachindex(fv)
        k = integ.problem.field_names[f]; v = copy(fv[f])
        v[c] = fields !== nothing && haskey(fields, k) ? _values_from(getproperty(fields, k), X) :
               _resample(v[c], length(X.φ))
        fv[f] = v
    end
    _set_fields!(integ, fv)
    _invalidate_velocities!(integ)
end

"""
    add_curve!(integ, curve; tag = nothing, fields = nothing) -> id

Add a curve (a `ThetaLState` or a counterclockwise parametrization) at the end of the list; returns its new id. Its
fields: from `fields` (a named tuple: per field a number, a vector or a function `X -> values`), otherwise from each
field's `init` (a number or a function).
"""
function add_curve!(integ::ThetaLIntegrator, curve; tag=nothing, fields=nothing)
    id = integ.next_id; integ.next_id += 1
    X = _as_state(integ, curve)
    fv = _field_vectors(integ)
    for f in eachindex(fv)
        k = integ.problem.field_names[f]
        v = fields !== nothing && haskey(fields, k) ? _values_from(getproperty(fields, k), X) :
            _field_values(integ.problem.fields[f], X, nothing)
        fv[f] = [fv[f]; [v]]
    end
    _set_fields!(integ, fv)
    integ.states = [integ.states; X]
    integ.ids = [integ.ids; id]; integ.tags = Any[integ.tags; tag]; integ.rigid = [integ.rigid; NO_RIGID_MOTION]
    _invalidate_velocities!(integ)
    id
end

"""Remove curve `c` (with its id and tag)."""
function remove_curve!(integ::ThetaLIntegrator, c::Int)
    keep = [i for i in eachindex(integ.states) if i != c]
    integ.states = integ.states[keep]; integ.ids = integ.ids[keep]; integ.tags = integ.tags[keep]
    integ.rigid = integ.rigid[keep]
    _set_fields!(integ, Vector{Vector{Float64}}[v[keep] for v in _field_vectors(integ)])
    _invalidate_velocities!(integ)
end

"""Move curve `c` to its parallel curve at distance `δ` along the outward normal (`offset_curve`)."""
offset_curve!(integ::ThetaLIntegrator, c::Int, δ::Real) = set_curve!(integ, c, offset_curve(integ.states[c], δ))

_as_state(integ::ThetaLIntegrator, curve) = curve isa ThetaLState ? curve : _from_curve(curve, integ.tol, 64, integ.Nmax)
