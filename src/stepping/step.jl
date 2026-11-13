# One step of every curve together, for any one-step scheme (src/stepping/schemes.jl): the stateless level under
# the integrator. Per curve and stage, the normal velocity is split into its rigid part (integrated exactly: it
# turns θ by Ω and moves x₀ on the rigid motion) and the deformation (Hou–Lowengrub–Shelley kinematics); φ, L and x₀
# are then advanced with the scheme's weights. One velocity evaluation per stage after the first, for all curves at
# once, and one at the end.
#
# Why: under a rigid motion φ̂_k (k ≠ 0) and L stay fixed and only the mean of φ moves (by Ωh); x₀ moves on the
# rigid motion, which `_rigid_flow` integrates exactly; and the exponential weights keep a state at rest fixed with
# a symbol present (e^{z}φ̂ − hϕ₁(z)σφ̂ = φ̂ with z = σh). So rigid translation and rotation are exact at any Δt.

"""
    StepResult

The outcome of `thetaL_step`: `status` (`:success` or `:rejected`), `reason` (why a step was rejected:
`:nonpositive_length`, `:error_estimate`, `:invalid_state`, or the reason of a `StepFailure`), the new `states`, their
velocities `U`, the context `ctx` of the last evaluation, the `rigid` motion of each curve over the step, the
embedded `error_estimate` (`NaN` if the scheme has none), the number of velocity `evaluations` made, and the state
carried with the curves: `fields` (one vector per field and curve, `fields[f][c]`), `bulk` and the last evaluation's
W, rates and ydot (`aux`). On rejection the states, velocities, context and carried state are those given.
"""
struct StepResult
    status::Symbol
    reason::Any
    states::Vector{ThetaLState}
    U::Vector{Vector{Float64}}
    ctx::Any
    rigid::Vector{RigidMotion}
    error_estimate::Float64
    evaluations::Int
    fields::Vector{Vector{Vector{Float64}}}
    bulk::Union{Nothing,Vector{Float64}}
    aux::_Aux
end

# What one stage contributes to the step, for one curve: the explicit remainder R = θ_t − σφ̂ (half spectrum), L_t, and
# the rigid motion as the coefficients of z′ = iΩz + a (complex x₀), and the deformation's drift d of x₀ with the rigid
# rotation Θ the stage's state has turned through since the start of the step (d turns with the curve: e^{−iΘ}d is
# the drift in the frame of the start).
struct _StageRecord
    R::Vector{ComplexF64}
    Lt::Float64
    motion::RigidMotion
    Ω::Float64
    a::ComplexF64
    d::ComplexF64
    Θ::Float64
end

# (φ̂ from the state's cache: no FFT of φ). With `T` given, the frame's tangential velocity is written into it.
function _stage_record(X::ThetaLState, U::Vector{Float64}, σ::Vector{Float64}, filt::Vector{Float64}, rigid::Symbol,
                       T::Union{Nothing,Vector{Float64}}=nothing, Θ::Float64=0.0,
                       W::Union{Nothing,Vector{Float64}}=nothing)
    w = _work(length(X.φ)); Ud = _take_r(w)
    motion = _rigid_part!(Ud, X, U, rigid, W)            # (from u = U n + W τ when the velocity gives W)
    R = _take_h(w)                                      # (handed back at the end of the step)
    Lt, drift = _thetaL_rhs_hat!(R, X, Ud, filt, T)
    T === nothing || _frame_tangential!(T, X, motion)
    Ω = motion.rotation_rate
    R[1] += w.N * Ω                                     # the rigid rotation turns every tangent angle (k = 0)
    φh = _phihat(X)
    for j in eachindex(R)
        R[j] -= σ[j] * φh[j]
    end
    _give!(w, Ud)
    V, xc = motion.translation, motion.center
    a = complex(V[1], V[2]) - im * Ω * complex(xc[1], xc[2])
    _StageRecord(R, Lt, motion, Ω, a, complex(drift[1], drift[2]), Θ)
end

# the material tangential velocity of curve c, if the evaluation gave one (a problem with fields transported along the
# curves requires it)
_curve_W(aux::_Aux, c::Int) = aux.W === nothing ? nothing : aux.W[c]

# x₀ after a time h on z′ = iΩz + a (exact for constant Ω, a): z ← e^{iΩh} z + h ϕ₁(iΩh) a
function _rigid_flow(x0::Vec2, Ω::Float64, a::ComplexF64, h::Float64)
    z = complex(x0[1], x0[2])
    w = cis(Ω * h) * z + h * _ϕ12(complex(0.0, Ω * h))[1] * a
    Vec2(real(w), imag(w))
end

# The weights of a step for one curve, per mode, as vectors over the half spectrum (z = σ_k h is fixed for the
# step, so they are computed once): E[i] = e^{cᵢz} (i ≤ ns the stages, ns + 1 the end), A[i][k], b[k], b̂[k].
struct _StepWeights
    E::Vector{Vector{Float64}}
    A::Vector{Vector{Vector{Float64}}}
    b::Vector{Vector{Float64}}
    b̂::Vector{Vector{Float64}}
end
# (vectors from the work arrays `w` when given, handed back by `_give_weights!`)
function _step_weights(scheme::OneStepScheme, σ::Vector{Float64}, h::Float64, w::Union{Nothing,_GridWork}=nothing)
    nodes = stage_nodes(scheme); ns = length(nodes); n = length(σ)
    E = Vector{Vector{Float64}}(undef, ns + 1); A = Vector{Vector{Vector{Float64}}}(undef, ns)
    b = Vector{Vector{Float64}}(undef, ns); b̂ = Vector{Vector{Float64}}(undef, ns)
    for i in 1:ns+1
        E[i] = _weight_vector(w, n)
    end
    for i in 1:ns
        A[i] = Vector{Vector{Float64}}(undef, i - 1)
        for k in 1:i-1
            A[i][k] = _weight_vector(w, n)
        end
        b[i] = _weight_vector(w, n); b̂[i] = _weight_vector(w, n)
    end
    for j in 1:n
        z = σ[j] * h
        wA, wb, wb̂ = scheme_weights(scheme, z)
        _store_rows!(A, j, wA, 1); _store!(b, j, wb); _store!(b̂, j, wb̂)
        for i in 1:ns
            E[i][j] = exp(nodes[i] * z)
        end
        E[ns+1][j] = exp(z)
    end
    _StepWeights(E, A, b, b̂)
end
_weight_vector(w::Nothing, n::Int) = Vector{Float64}(undef, n)
_weight_vector(w::_GridWork, n::Int) = _take_hr(w)
function _give_weights!(w::_GridWork, W::_StepWeights)
    for v in W.E
        _give!(w, v)
    end
    for i in eachindex(W.b)
        _give!(w, W.b[i]); _give!(w, W.b̂[i])
        for v in W.A[i]
            _give!(w, v)
        end
    end
end
# entry j of the vectors `dest[k]` from the tuple t (homogeneous); the rows of A by recursion on the tuple's type
# (rows of different lengths: a loop over them would not be type-stable)
function _store!(dest::Vector{Vector{Float64}}, j::Int, t::Tuple)
    for k in 1:length(t)
        dest[k][j] = t[k]
    end
end
_store_rows!(A::Vector{Vector{Vector{Float64}}}, j::Int, t::Tuple{}, i::Int) = nothing
_store_rows!(A::Vector{Vector{Vector{Float64}}}, j::Int, t::Tuple, i::Int) =
    (_store!(A[i], j, first(t)); _store_rows!(A, j, Base.tail(t), i + 1))

# The state reached from `X` (half spectrum φh) with stage records `rec[1:m]`, over the time c·h: the weights `E`
# (e^{cz} per mode) and `W[k]` (per record, per mode), and `W0[k]` the same at z = 0, for L and x₀ (they sum to c:
# the rigid motion averaged over the stages). Writes the filtered φ̂ into `out` (a half spectrum; the coefficients
# of k = 0 and N/2 made real, as those of a real function are) and returns it, L and x₀. `c = 1` is the end.
function _combine!(out::Vector{ComplexF64}, X::ThetaLState, φh::Vector{ComplexF64}, E::Vector{Float64},
                   W::Vector{Vector{Float64}}, W0::Vector{Vector{Float64}}, filt::Vector{Float64},
                   rec::Vector{_StageRecord}, m::Int, c::Float64, h::Float64)
    n = length(φh); N = length(X.φ)
    for j in 1:n
        acc = 0.0im
        for k in 1:m
            acc += W[k][j] * rec[k].R[j]
        end
        out[j] = (E[j] * φh[j] + h * acc) * filt[j]
    end
    out[1] = real(out[1]); iseven(N) && (out[n] = real(out[n]))
    Lacc = 0.0; Ω = 0.0; a = 0.0im; D = 0.0im
    for k in 1:m
        w = W0[k][1]
        Lacc += w * rec[k].Lt; Ω += (w / c) * rec[k].Ω; a += (w / c) * rec[k].a
        D += w * cis(-rec[k].Θ) * rec[k].d                 # the drifts in the frame of the start
    end
    # x₀ on the rigid motion, plus the deformation's drift turned with the curve: z′ = iΩz + a + e^{iΩt}D/c
    x0 = _rigid_flow(X.x0, Ω, a, c * h); Θ = Ω * c * h; δ = h * cis(Θ) * D
    out, X.L + h * Lacc, x0 + Vec2(real(δ), imag(δ)), Θ
end

# The state with half spectrum φh, length L and point x0, φh kept as its cached spectrum
function _state_from_spectrum(φh::Vector{ComplexF64}, N::Int, L::Float64, x0::Vec2)
    w = _work(N); c = copyto!(_take_h(w), φh)
    φ = _irfft!(w, Vector{Float64}(undef, N), c); _give!(w, c)
    ThetaLState(φ, L, x0, _StateCache(φh, nothing, nothing))
end

# The same, with the mean of φ brought back to (−π, π] (a full turn changes nothing; the cached spectrum too)
function _wrapped_state(φh::Vector{ComplexF64}, N::Int, L::Float64, x0::Vec2)
    X = _state_from_spectrum(φh, N, L, x0); φ = X.φ
    n = round(sum(φ) / (2π * N))
    if n != 0                                    # (X is not yet shared: changing it here is safe)
        for j in 1:N
            φ[j] -= 2π * n
        end
        φh[1] -= 2π * n * N
    end
    X
end

"""
    thetaL_step(scheme, problem, states, U, t, t_end; warm = nothing, ids, tags, topology = 0, attempt = 1,
                errtol = Inf, fields = nothing, bulk = nothing, aux = nothing) -> StepResult

One step of all `states` from `t` to `t_end` (exactly) with the scheme, the velocity, symbol, rigid mode and
`validate` of `problem`, starting from the velocities `U` at `t` and the context `warm`. With state carried with the
curves: its values at `t` (`fields[f][c]`, `bulk`; by default the problem's initial values) and the W, rates and
ydot of the evaluation that gave `U` (`aux`, as a velocity's named tuple; required then). Stateless: the integrator
(`advance!`) is built on it. A step is rejected (nothing changed) when a stage has a non-positive length, when the
embedded error estimate exceeds `errtol` (checked before the end evaluation), when `validate` fails, or when the
velocity throws `StepFailure`.
"""
function thetaL_step(scheme::OneStepScheme, problem::CurveProblem, states::Vector{ThetaLState},
                     U::Vector{Vector{Float64}}, t::Real, t_end::Real; warm=nothing,
                     ids::Vector{Int}=collect(1:length(states)), tags::Vector{Any}=Any[nothing for _ in states],
                     topology::Int=0, attempt::Int=1, errtol::Real=Inf, fields=nothing, bulk=nothing, aux=nothing)
    fv = fields === nothing ? _initial_fields(problem, states) : fields
    by = bulk === nothing && problem.bulk !== nothing ? copy(problem.bulk.y0) : bulk
    ax = aux === nothing ? _NO_AUX : aux isa _Aux ? aux : _read_aux(problem, aux, states)
    _thetaL_step(scheme, problem, states, U, Float64(t), Float64(t_end), warm, ids, tags, topology, attempt,
                 Float64(errtol), fv, by, ax)
end

_initial_fields_tuple(problem::CurveProblem, states::Vector{ThetaLState}) =
    _fields_tuple(problem.field_names, _initial_fields(problem, states))

# the problem's initial field values on the states: fields[f][c]
function _initial_fields(problem::CurveProblem, states::Vector{ThetaLState})
    fv = Vector{Vector{Vector{Float64}}}(undef, length(problem.fields))
    for f in eachindex(problem.fields)
        fv[f] = Vector{Vector{Float64}}(undef, length(states))
        for c in eachindex(states)
            fv[f][c] = _field_values(problem.fields[f], states[c], c)
        end
    end
    fv
end

# A function called through this is not inferred from the caller: the carried state's code is compiled only for
# problems that have it (precompile time), at the cost of one dynamic call per stage
_barrier(f) = Base.inferencebarrier(f)
const _NO_FIELDS = NamedTuple()

Base.@nospecializeinfer function _thetaL_step(scheme::OneStepScheme, problem::CurveProblem,
                                              states::Vector{ThetaLState}, U::Vector{Vector{Float64}},
                                              t::Float64, t_end::Float64, warm, ids::Vector{Int}, tags::Vector{Any},
                                              topology::Int, attempt::Int, errtol::Float64,
                                              fields::Vector{Vector{Vector{Float64}}},
                                              bulk::Union{Nothing,Vector{Float64}}, aux::_Aux)
    @nospecialize warm
    h = t_end - t; n = length(states); nodes = stage_nodes(scheme); ns = length(nodes)
    reject(reason, err, evals) = StepResult(:rejected, reason, states, U, warm, RigidMotion[], err, evals, fields,
                                            bulk, aux)
    carried = !isempty(problem.fields) || problem.bulk !== nothing
    cs = carried ? _barrier(_carried_start)(scheme, problem, states, fields, bulk, aux, h) : nothing
    T = carried ? cs.T::Union{Nothing,Vector{Vector{Float64}}} : nothing   # the frame's T, if a field needs it

    # the start: filter, frozen symbol and the step's weights, Fourier φ and the first stage's record, per curve
    filt = Vector{Vector{Float64}}(undef, n); σ = similar(filt); φh = Vector{Vector{ComplexF64}}(undef, n)
    W = Vector{_StepWeights}(undef, n); W0 = _step_weights(scheme, [0.0], h)      # (W0: the weights at z = 0)
    rec = [Vector{_StageRecord}(undef, ns) for _ in 1:n]; Θ = [zeros(ns) for _ in 1:n]   # (Θ: each stage's rotation)
    for c in 1:n
        N = length(states[c].φ)
        filt[c] = _cached_filter(N, scheme.filter)
        σ[c] = _curve_rates(problem.symbol, states[c], _work(N).absk)
        W[c] = _step_weights(scheme, σ[c], h, _work(N))
        φh[c] = _phihat(states[c])
        rec[c][1] = _stage_record(states[c], U[c], σ[c], filt[c], problem.rigid, T === nothing ? nothing : T[c], 0.0,
                                  _curve_W(aux, c))
    end
    carried && _barrier(_carried_records!)(cs, problem, 1, states, rec, aux)

    # the stages after the first: one evaluation each, all curves together
    ctx = warm; evaluations = 0
    for i in 2:ns
        stage_states = Vector{ThetaLState}(undef, n)
        for c in 1:n
            N = length(states[c].φ)
            φi, Li, x0i, Θ[c][i] = _combine!(Vector{ComplexF64}(undef, N ÷ 2 + 1), states[c], φh[c], W[c].E[i], W[c].A[i],
                                    W0.A[i], filt[c], rec[c], i - 1, nodes[i], h)
            Li > 0 || return reject(:nonpositive_length, NaN, evaluations)
            stage_states[c] = _state_from_spectrum(φi, N, Li, x0i)
        end
        fi, yi = carried ? _barrier(_carried_stage!)(cs, i, h, filt, stage_states) : (_NO_FIELDS, nothing)
        ev = Evaluation(t + nodes[i] * h, :stage, i, nodes[i], h, attempt, ids, tags, topology, states,
                        fi::NamedTuple, yi::Union{Nothing,Vector{Float64}})
        Ui, auxi = try
            Ui, ctx, auxi = _evaluate_velocity(problem, stage_states, ev, ctx); Ui, auxi
        catch e
            e isa StepFailure || rethrow()
            return reject(e.reason, NaN, evaluations + 1)
        end
        evaluations += 1
        for c in 1:n
            rec[c][i] = _stage_record(stage_states[c], Ui[c], σ[c], filt[c], problem.rigid,
                                      T === nothing ? nothing : T[c], Θ[c][i], _curve_W(auxi, c))
        end
        carried && _barrier(_carried_records!)(cs, problem, i, stage_states, rec, auxi)
    end

    # the end state, and the embedded error estimate before the end evaluation
    end_states = Vector{ThetaLState}(undef, n); err = has_error_estimate(scheme) ? 0.0 : NaN
    for c in 1:n
        N = length(states[c].φ); w = _work(N)
        φn, Ln, x0n, _ = _combine!(Vector{ComplexF64}(undef, w.n), states[c], φh[c], W[c].E[ns+1], W[c].b, W0.b, filt[c],
                                rec[c], ns, 1.0, h)
        Ln > 0 || return reject(:nonpositive_length, NaN, evaluations)
        if has_error_estimate(scheme)
            φe, Le, x0e, _ = _combine!(_take_h(w), states[c], φh[c], W[c].E[ns+1], W[c].b̂, W0.b̂, filt[c], rec[c], ns,
                                    1.0, h)
            for j in eachindex(φe)
                φe[j] = φn[j] - φe[j]
            end
            dφ = _irfft!(w, _take_r(w), φe)
            err = max(err, maximum(abs, dφ), abs(Ln - Le) / Ln, sqrt(sum(abs2, x0n - x0e)) / Ln)
            _give!(w, φe, dφ)
        end
        end_states[c] = _wrapped_state(φn, N, Ln, x0n)
    end
    fn = fields; yn = bulk
    if carried
        fn_, yn_, errc = _barrier(_carried_end)(cs, scheme, h, filt, end_states)
        fn = fn_::Vector{Vector{Vector{Float64}}}; yn = yn_::Union{Nothing,Vector{Float64}}
        has_error_estimate(scheme) && (err = max(err, errc::Float64))
    end
    err > errtol && return reject(:error_estimate, err, evaluations)
    problem.validate === nothing || problem.validate(end_states) || return reject(:invalid_state, err, evaluations)

    # the end evaluation (also the next step's first stage)
    fnt = carried ? _barrier(_fields_tuple)(problem.field_names, fn)::NamedTuple : _NO_FIELDS
    ev = Evaluation(t_end, :end, ns + 1, 1.0, h, attempt, ids, tags, topology, states, fnt, yn)
    Un, auxn = try
        Un, ctx, auxn = _evaluate_velocity(problem, end_states, ev, ctx); Un, auxn
    catch e
        e isa StepFailure || rethrow()
        return reject(e.reason, err, evaluations + 1)
    end
    evaluations += 1

    # the rigid motion over the step: the stages' motions with the end weights at z = 0
    rigid = Vector{RigidMotion}(undef, n)
    for c in 1:n
        V = Vec2(0.0, 0.0); Ω = 0.0; xc = Vec2(0.0, 0.0)
        for k in 1:ns
            m = rec[c][k].motion; w0 = W0.b[k][1]
            V += w0 * m.translation; Ω += w0 * m.rotation_rate; xc += w0 * m.center
        end
        rigid[c] = RigidMotion(V, Ω, xc)
    end
    for c in 1:n                                        # the step's work vectors back to the work arrays
        w = _work(length(states[c].φ)); _give_weights!(w, W[c])
        for r in rec[c]
            _give!(w, r.R)
        end
    end
    carried && _barrier(_carried_release!)(cs, states)
    StepResult(:success, nothing, end_states, Un, ctx, rigid, err, evaluations, fn, yn, auxn)
end
