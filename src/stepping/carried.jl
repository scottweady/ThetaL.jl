# The step's part for state carried with the curves (problem/carried.jl): per stage, the explicit remainder of each
# field on each curve and of the bulk state, and their combination with the scheme's weights, as for φ.
#
# A field c on a curve is advanced on its half spectrum like φ: u_t ≈ σ û for its own symbol σ (frozen per curve
# and step), the rest R = u_t − σû explicit, filtered as φ is. The velocity's rate S is completed by the transport
# along the curve: the material moves tangentially at W, the grid points at the frame's T (the rigid motion's
# tangential velocity plus the equal-arclength one of the deformation, T_d), so the field slides along the grid at
# W − T. The variable u advanced:
#   :density   u = m = s_α c, the density per unit α:  m_t = s_α S − ((W − T)c)_α
#              (c_t = S − ((W − T)c)_s − c L_t/L; with m the L_t term cancels and Σm is a linear invariant, which
#              every Runge–Kutta step keeps to roundoff: ∫c ds = (2π/N) Σm exactly; c = m/s_α at each stage)
#   :material  u = c:  c_t = S − (W − T) c_α/s_α
#   :none      u = c:  c_t = S
# The bulk state y is advanced the same way in its own components: y_t ≈ σ_y y, R = ydot − σ_y y.

# The tangential velocity of the grid points of X: the rigid motion's (translation V and rotation Ω about x_c:
# V·τ + Ω (x − x_c)·n) plus T_d (in T on entry), written into T
function _frame_tangential!(T::Vector{Float64}, X::ThetaLState, motion::RigidMotion)
    N = length(X.φ); x = _points(X); V = motion.translation; Ω = motion.rotation_rate; xc = motion.center
    for j in 1:N
        sn, cs = sincos(_α(j, N) + X.φ[j]); r = x[j] - xc
        T[j] += V[1] * cs + V[2] * sn + Ω * (r[1] * sn - r[2] * cs)
    end
    T
end

# The explicit remainder of a field at a stage (half spectrum, from the work arrays): values cv on curve X and the
# spectrum ch of the variable advanced (m = s_α c for :density, c otherwise), the field's `transport`, the material
# and frame tangential velocities Wt and T (unused for :none), the velocity's rate S (`nothing`: 0) and the field's
# frozen symbol σ
function _field_record(X::ThetaLState, cv::Vector{Float64}, ch::Vector{ComplexF64}, transport::Symbol,
                       Wt::Union{Nothing,Vector{Float64}}, T::Union{Nothing,Vector{Float64}},
                       S::Union{Nothing,Vector{Float64}}, σ::Vector{Float64})
    N = length(X.φ); w = _work(N); sα = X.L / 2π
    g = _take_r(w); R = _take_h(w)
    for j in 1:N
        g[j] = S === nothing ? 0.0 : S[j]
    end
    if transport === :density                           # m_t = s_α S − ((W − T)c)_α
        flux = _take_r(w)
        for j in 1:N
            flux[j] = (Wt[j] - T[j]) * cv[j]; g[j] *= sα
        end
        _rfft!(w, R, g); F = _rfft!(w, _take_h(w), flux)
        for j in eachindex(R)
            k = j - 1
            _nyquist(k, N) || (R[j] -= im * k * F[j])
        end
        _give!(w, flux, F)
    elseif transport === :material
        cα = _irfft!(w, _take_r(w), _dα_half!(copyto!(_take_h(w), ch), N))
        for j in 1:N
            g[j] -= (Wt[j] - T[j]) * cα[j] / sα
        end
        _rfft!(w, R, g); _give!(w, cα)
    else
        _rfft!(w, R, g)
    end
    for j in eachindex(R)
        R[j] -= σ[j] * ch[j]
    end
    _give!(w, g)
    R
end

# A field's half spectrum after the time c·h (as `_combine!` for φ): out = (E ĉ + h Σ_k W[k] R[k]) · filt, the
# coefficients of k = 0 and N/2 made real
function _combine_field!(out::Vector{ComplexF64}, ch::Vector{ComplexF64}, E::Vector{Float64},
                         W::Vector{Vector{Float64}}, R::Vector{Vector{ComplexF64}}, m::Int, h::Float64,
                         filt::Vector{Float64}, N::Int)
    n = length(ch)
    for j in 1:n
        acc = 0.0im
        for k in 1:m
            acc += W[k][j] * R[k][j]
        end
        out[j] = (E[j] * ch[j] + h * acc) * filt[j]
    end
    out[1] = real(out[1]); iseven(N) && (out[n] = real(out[n]))
    out
end

# The bulk state after the time c·h: out = E y + h Σ_k W[k] R[k], componentwise
function _combine_bulk!(out::Vector{Float64}, y::Vector{Float64}, E::Vector{Float64}, W::Vector{Vector{Float64}},
                        R::Vector{Vector{Float64}}, m::Int, h::Float64)
    for j in eachindex(y)
        acc = 0.0
        for k in 1:m
            acc += W[k][j] * R[k][j]
        end
        out[j] = E[j] * y[j] + h * acc
    end
    out
end

# the bulk remainder ydot − σ_y y
function _bulk_record(ydot::Vector{Float64}, σy::Vector{Float64}, y::Vector{Float64})
    R = similar(y)
    for j in eachindex(y)
        R[j] = ydot[j] - σy[j] * y[j]
    end
    R
end

# the largest |a − b| relative to the largest |a| (0 for a = b = 0)
function _rel_diff(a::Vector{Float64}, b::Vector{Float64})
    d = 0.0; s = 0.0
    for j in eachindex(a)
        d = max(d, abs(a[j] - b[j])); s = max(s, abs(a[j]))
    end
    d == 0 ? 0.0 : d / s
end

# The carried state through one step, kept apart from the curves' step: `_thetaL_step` calls these through an
# inference barrier, so a problem without carried state never compiles them (precompile time) and one with it
# compiles them at first use. Per field f and curve c: the frozen symbol, the step's weights, the spectrum at the
# start and the stage records; for the bulk the same; and the frame's tangential velocity T per curve, written by
# the curves' stage records when a field is transported.
struct _CarriedStep
    names::Vector{Symbol}
    density::Vector{Bool}                                 # per field: advanced as m = s_α c
    fields::Vector{Vector{Vector{Float64}}}               # the values at the start, fields[f][c]
    bulk::Union{Nothing,Vector{Float64}}
    σf::Vector{Vector{Vector{Float64}}}
    Wf::Vector{Vector{_StepWeights}}
    chf::Vector{Vector{Vector{ComplexF64}}}
    recf::Vector{Vector{Vector{Vector{ComplexF64}}}}      # recf[f][c][stage]
    σy::Vector{Float64}
    Wy::Union{Nothing,_StepWeights}
    recy::Vector{Vector{Float64}}
    T::Union{Nothing,Vector{Vector{Float64}}}
    stage_fields::Vector{Vector{Vector{Float64}}}        # the current stage's values and spectra
    stage_ch::Vector{Vector{Vector{ComplexF64}}}
    stage_bulk::Vector{Union{Nothing,Vector{Float64}}}
end

function _carried_start(scheme::OneStepScheme, problem::CurveProblem, states::Vector{ThetaLState},
                        fields::Vector{Vector{Vector{Float64}}}, bulk::Union{Nothing,Vector{Float64}}, aux::_Aux,
                        h::Float64)
    n = length(states); nf = length(problem.fields); ns = length(stage_nodes(scheme))
    transported = any(f -> f.transport !== :none, problem.fields)
    length(fields) == nf && all(length(fv) == n for fv in fields) ||
        throw(DimensionMismatch("the fields' values do not match the problem's fields and the curves"))
    nf > 0 && transported && aux.W === nothing &&
        throw(ArgumentError("a field is transported along the curves: the step needs W (the velocity returns it)"))
    problem.bulk !== nothing && (bulk === nothing || aux.ydot === nothing) &&
        throw(ArgumentError("the problem has a bulk state: the step needs its value and ydot"))
    σf = [Vector{Vector{Float64}}(undef, n) for _ in 1:nf]; Wf = [Vector{_StepWeights}(undef, n) for _ in 1:nf]
    chf = [Vector{Vector{ComplexF64}}(undef, n) for _ in 1:nf]
    recf = [[Vector{Vector{ComplexF64}}(undef, ns) for _ in 1:n] for _ in 1:nf]
    density = [f.transport === :density for f in problem.fields]
    for f in 1:nf, c in 1:n
        N = length(states[c].φ)
        σf[f][c] = _curve_rates(problem.fields[f].symbol, states[c], _work(N).absk)
        Wf[f][c] = _step_weights(scheme, σf[f][c], h, _work(N))
        chf[f][c] = _rfft(density[f] ? fields[f][c] .* (states[c].L / 2π) : fields[f][c])   # m = s_α c, or c
    end
    hasbulk = problem.bulk !== nothing
    σy = hasbulk ? (problem.bulk.symbol === nothing ? zeros(length(bulk)) : problem.bulk.symbol) : Float64[]
    Wy = hasbulk ? _step_weights(scheme, σy, h) : nothing
    T = transported && nf > 0 ? [Vector{Float64}(undef, length(X.φ)) for X in states] : nothing
    _CarriedStep(copy(problem.field_names), density, fields, bulk, σf, Wf, chf, recf, σy, Wy,
                 Vector{Vector{Float64}}(undef, ns), T, copy(fields), copy(chf), Union{Nothing,Vector{Float64}}[bulk])
end

# the records of stage i from the stage's states, the curves' records (for L_t) and the evaluation's aux
function _carried_records!(cs::_CarriedStep, problem::CurveProblem, i::Int, states::Vector{ThetaLState},
                           rec::Vector{Vector{_StageRecord}}, aux::_Aux)
    for f in eachindex(cs.names), c in eachindex(states)
        cs.recf[f][c][i] = _field_record(states[c], cs.stage_fields[f][c], cs.stage_ch[f][c],
                                         problem.fields[f].transport, aux.W === nothing ? nothing : aux.W[c],
                                         cs.T === nothing ? nothing : cs.T[c],
                                         aux.rates[f] === nothing ? nothing : aux.rates[f][c], cs.σf[f][c])
    end
    cs.Wy === nothing || (cs.recy[i] = _bulk_record(aux.ydot, cs.σy, cs.stage_bulk[1]))
    cs
end

# the values at stage i (from the records of the stages before it) on the stage's curves: the fields as the
# velocity sees them (c = m/s_α for a density, with the stage's length), the bulk
function _carried_stage!(cs::_CarriedStep, i::Int, h::Float64, filt::Vector{Vector{Float64}},
                         states::Vector{ThetaLState})
    n = length(states); nf = length(cs.names)
    fi = [Vector{Vector{Float64}}(undef, n) for _ in 1:nf]
    chi = [Vector{Vector{ComplexF64}}(undef, n) for _ in 1:nf]
    for f in 1:nf, c in 1:n
        N = length(states[c].φ); W = cs.Wf[f][c]
        chi[f][c] = _combine_field!(Vector{ComplexF64}(undef, N ÷ 2 + 1), cs.chf[f][c], W.E[i], W.A[i],
                                    cs.recf[f][c], i - 1, h, filt[c], N)
        fi[f][c] = _irfft(chi[f][c], N)
        cs.density[f] && (fi[f][c] .*= 2π / states[c].L)
    end
    yi = cs.Wy === nothing ? nothing :
         _combine_bulk!(similar(cs.bulk), cs.bulk, cs.Wy.E[i], cs.Wy.A[i], cs.recy, i - 1, h)
    cs.stage_fields .= fi; cs.stage_ch .= chi; cs.stage_bulk[1] = yi
    _fields_tuple(cs.names, fi), yi
end

# the values at the end of the step on the end curves and their part of the embedded error estimate (relative to
# their size)
function _carried_end(cs::_CarriedStep, scheme::OneStepScheme, h::Float64, filt::Vector{Vector{Float64}},
                      states::Vector{ThetaLState})
    n = length(states); nf = length(cs.names); ns = length(stage_nodes(scheme)); err = 0.0
    est = has_error_estimate(scheme)
    fn = [Vector{Vector{Float64}}(undef, n) for _ in 1:nf]
    for f in 1:nf, c in 1:n
        N = length(states[c].φ); w = _work(N); W = cs.Wf[f][c]
        cn = _combine_field!(_take_h(w), cs.chf[f][c], W.E[ns+1], W.b, cs.recf[f][c], ns, h, filt[c], N)
        fn[f][c] = _irfft(cn, N); s = cs.density[f] ? 2π / states[c].L : 1.0
        s == 1.0 || (fn[f][c] .*= s)
        if est
            ce = _combine_field!(_take_h(w), cs.chf[f][c], W.E[ns+1], W.b̂, cs.recf[f][c], ns, h, filt[c], N)
            de = _irfft!(w, _take_r(w), ce); s == 1.0 || (de .*= s)
            err = max(err, _rel_diff(fn[f][c], de))
            _give!(w, ce, de)
        end
        _give!(w, cn)
    end
    yn = nothing
    if cs.Wy !== nothing
        yn = _combine_bulk!(similar(cs.bulk), cs.bulk, cs.Wy.E[ns+1], cs.Wy.b, cs.recy, ns, h)
        est && (err = max(err, _rel_diff(yn, _combine_bulk!(similar(cs.bulk), cs.bulk, cs.Wy.E[ns+1], cs.Wy.b̂,
                                                             cs.recy, ns, h))))
    end
    fn, yn, err
end

# the step's work vectors back to the work arrays
function _carried_release!(cs::_CarriedStep, states::Vector{ThetaLState})
    for f in eachindex(cs.names), c in eachindex(states)
        w = _work(length(states[c].φ)); _give_weights!(w, cs.Wf[f][c])
        for R in cs.recf[f][c]
            _give!(w, R)
        end
    end
    cs
end
