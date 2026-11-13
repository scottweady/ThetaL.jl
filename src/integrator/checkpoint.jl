# Checkpoints: the state of a run as plain data (numbers, vectors, strings and the caller's tags), so that it can be
# written with any serializer (JLD2, HDF5, Serialization) and read back by a later version of ThetaL. Not saved:
# the problem's functions and the velocity context (a restart evaluates the velocities again only if U is absent).

const CHECKPOINT_FORMAT = 1

"""
    checkpoint(integ) -> NamedTuple

The state of the run as plain data: `format` (a version number), `t`, per curve `phi`, `phi_spectrum` (the half
spectrum of φ the step works with, so that a restart continues bit for bit), `length`, `x0` (as `[x, y]`), `U` (or
`nothing`), `ids`, `tags`, `topology`, `next_id`, the `scheme` (`type`, `filter_order`, `filter_cutoff`, resolved),
the `stats` counts, and the state carried with the curves: `fields` (a named tuple, one vector per curve each),
`bulk`, and with U the evaluation's `W`, `rates` and `ydot`. Restart with `init(problem, scheme; from = checkpoint)`; a checkpoint without
`phi_spectrum` is read too, its spectra computed from φ.
"""
function checkpoint(integ::ThetaLIntegrator)
    f = integ.scheme.filter
    (; format=CHECKPOINT_FORMAT, t=integ.t,
       phi=[copy(X.φ) for X in integ.states], phi_spectrum=[copy(_phihat(X)) for X in integ.states],
       length=[X.L for X in integ.states],
       x0=[[X.x0[1], X.x0[2]] for X in integ.states],
       U=integ.U === nothing ? nothing : [copy(u) for u in integ.U],
       ids=copy(integ.ids), tags=copy(integ.tags), topology=integ.topology, next_id=integ.next_id,
       scheme=(; type=string(nameof(typeof(integ.scheme))), filter_order=f.order, filter_cutoff=f.cutoff),
       stats=(; steps=integ.stats.steps, evaluations=integ.stats.evaluations, rejections=integ.stats.rejections),
       fields=_fields_tuple(integ.problem.field_names, [[copy(v) for v in fv] for fv in _field_vectors(integ)]),
       bulk=integ.bulk === nothing ? nothing : copy(integ.bulk),
       W=integ.aux.W === nothing ? nothing : [copy(v) for v in integ.aux.W],
       rates=_fields_tuple(integ.problem.field_names,
                           [r === nothing ? nothing : [copy(v) for v in r] for r in integ.aux.rates]),
       ydot=integ.aux.ydot === nothing ? nothing : copy(integ.aux.ydot))
end

function _init_from_checkpoint(problem::CurveProblem, scheme::OneStepScheme, cp, tol::Float64, Nmax::Int)
    cp.format == CHECKPOINT_FORMAT || throw(ArgumentError("checkpoint format $(cp.format); this ThetaL reads $CHECKPOINT_FORMAT"))
    states = ThetaLState[ThetaLState(Vector{Float64}(cp.phi[c]), Float64(cp.length[c]), Vec2(cp.x0[c][1], cp.x0[c][2]))
                         for c in eachindex(cp.phi)]
    if haskey(cp, :phi_spectrum)                     # the spectra the run was using: the restart continues exactly
        for c in eachindex(states)
            states[c].cache.φh = Vector{ComplexF64}(cp.phi_spectrum[c])
        end
    end
    f = scheme.filter                                # an :auto cutoff takes the one resolved in the first run
    if f.cutoff === :auto && cp.scheme.filter_cutoff isa Integer
        scheme = with_filter(scheme, SpectralFilter(order=f.order, cutoff=Int(cp.scheme.filter_cutoff)))
    end
    stats = IntegratorStats(cp.stats.steps, cp.stats.evaluations, cp.stats.rejections, 0, NaN, nothing, NaN)
    names = problem.field_names
    carried = !isempty(names) || problem.bulk !== nothing
    carried && !haskey(cp, :fields) && throw(ArgumentError("the checkpoint has no carried state for this problem"))
    fv = Vector{Vector{Float64}}[]
    for k in names
        haskey(cp.fields, k) || throw(ArgumentError("the checkpoint has no field $k"))
        push!(fv, [Vector{Float64}(v) for v in getproperty(cp.fields, k)])
    end
    bulk = problem.bulk === nothing ? nothing : Vector{Float64}(cp.bulk)
    integ = ThetaLIntegrator(problem, scheme, Float64(cp.t), states, nothing, nothing, Vector{Int}(cp.ids),
                             Vector{Any}(cp.tags), Int(cp.topology), Int(cp.next_id), fill(NO_RIGID_MOTION, length(states)),
                             stats, tol, Nmax, _fields_tuple(names, fv), bulk, _NO_AUX)
    if cp.U !== nothing && carried                   # the evaluation's W, rates and ydot, as it gave them
        W = cp.W === nothing ? nothing : _checked_curve_vectors(cp.W, states, "W")
        rates = Union{Nothing,Vector{Vector{Float64}}}[
            haskey(cp.rates, k) && getproperty(cp.rates, k) !== nothing ?
            _checked_curve_vectors(getproperty(cp.rates, k), states, "rate of $k") : nothing for k in names]
        integ.aux = _Aux(W, rates, cp.ydot === nothing ? nothing : Vector{Float64}(cp.ydot))
    end
    cp.U === nothing ? _evaluate_initial!(integ) : (integ.U = _checked_velocities(cp.U, states))
    _resolve_auto_cutoff!(integ)
end
