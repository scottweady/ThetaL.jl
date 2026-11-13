# State carried with the curves: fields on the curves (`CurveField`: a value per grid point of each curve, e.g. a
# concentration c) and a bulk state off them (`BulkState`: one vector y, e.g. the unknowns of a field with memory).
# Both are advanced by the step with the curves, with the scheme's weights, rolled back on rejection, kept through
# curve edits and checkpoints. The velocity sees their stage values (`ev.fields`, `ev.bulk`) and returns their rates
# (`(; U, W, rates, ydot)`); the step itself is in src/stepping/carried.jl.

"""
    CurveField(init; transport = :density, symbol = nothing)

A field on the curves: one value per grid point of each curve, advanced with them.

- `init`: the initial values, for every curve the same (a number, or a function `X -> values` of the curve's
  state), or one entry per curve (each a number, a vector of the curve's N values, or such a function). A function
  or a number is also used for curves added later (`add_curve!`).
- `transport`: how the field moves along the curve, where the material moves tangentially at W (returned by the
  velocity) and the grid points at the frame's T (which keeps them at equal arclength):
  - `:density`: a conserved density per unit length; the step adds c_t = −((W − T)c)_s − c L_t/L, which keeps
    ∫c ds exactly on the grid (the velocity's rate is then the source, e.g. adsorption or surface diffusion);
  - `:material`: a value carried by the material points: c_t = −(W − T) c_s;
  - `:none`: the velocity's rate is the whole c_t at the grid points.
- `symbol`: the stiff part of the field's own rate, as the curve's symbol: the growth rate `σ(q)` or `σ(q, X)` of a
  mode of wavenumber q along the curve (surface diffusion D c_ss: `q -> -D * q^2`), treated exactly by the step.
"""
struct CurveField
    init::Any
    transport::Symbol
    symbol::Any
end
function CurveField(init; transport::Symbol=:density, symbol=nothing)
    transport in (:density, :material, :none) ||
        throw(ArgumentError("transport must be :density, :material or :none, got :$transport"))
    CurveField(init, transport, symbol)
end

"""
    BulkState(y0; symbol = nothing)

A state off the curves, one vector y advanced with them (the unknowns of a bulk field with memory, a few scalars).
Its rate is the velocity's `ydot`. `symbol`: the growth rates σ of its components, a vector (y_t ≈ σ y for the stiff
part, diagonal: a pseudo-spectral field's −Dk²), treated exactly by the step; `nothing` for none.
"""
struct BulkState
    y0::Vector{Float64}
    symbol::Union{Nothing,Vector{Float64}}
end
function BulkState(y0::AbstractVector{<:Real}; symbol=nothing)
    σ = symbol === nothing ? nothing : Vector{Float64}(symbol)
    σ === nothing || length(σ) == length(y0) ||
        throw(DimensionMismatch("the bulk symbol has $(length(σ)) rates for $(length(y0)) components"))
    BulkState(Vector{Float64}(y0), σ)
end

# the initial values of a field on curve c (of states X); `i` the curve's position among the problem's curves (for
# per-curve values) or `nothing` for a curve added later
Base.@nospecializeinfer function _field_values(f::CurveField, X::ThetaLState, i::Union{Nothing,Int})
    spec = f.init
    if spec isa AbstractVector                          # one entry per curve
        i === nothing && throw(ArgumentError("a field given per curve needs initial values for an added curve " *
                                             "(`add_curve!(…; fields = (; c = …))`)"))
        spec = spec[i]
    end
    _values_from(spec, X)
end
Base.@nospecializeinfer function _values_from(@nospecialize(spec), X::ThetaLState)
    N = length(X.φ)
    v = spec isa Real ? fill(Float64(spec), N) : spec isa AbstractVector ? Vector{Float64}(spec) :
        Vector{Float64}(spec(X))
    length(v) == N || throw(DimensionMismatch("a field has $(length(v)) values for a curve of $N points"))
    v
end

# What a velocity evaluation returns besides U and ctx: the material tangential velocity W (one vector per curve;
# `nothing` if not given), the fields' rates (one entry per field: one vector per curve, or `nothing` for zero) and
# the bulk rate ydot.
struct _Aux
    W::Union{Nothing,Vector{Vector{Float64}}}
    rates::Vector{Union{Nothing,Vector{Vector{Float64}}}}
    ydot::Union{Nothing,Vector{Float64}}
end
const _NO_AUX = _Aux(nothing, Union{Nothing,Vector{Vector{Float64}}}[], nothing)

# the stage values handed to the velocity: the fields as a named tuple (one vector per curve each)
_fields_tuple(names::Vector{Symbol}, values::AbstractVector) = NamedTuple{Tuple(names)}(Tuple(values))

# per-curve vectors, checked against the curves' sizes, copied
Base.@nospecializeinfer function _checked_curve_vectors(@nospecialize(V), states::Vector{ThetaLState}, what::String)
    length(V) == length(states) ||
        throw(DimensionMismatch("the velocity returned $(length(V)) vectors of $what for $(length(states)) curves"))
    out = Vector{Vector{Float64}}(undef, length(states))
    for c in eachindex(states)
        v = Vector{Float64}(V[c])
        length(v) == length(states[c].φ) ||
            throw(DimensionMismatch("curve $c has $(length(states[c].φ)) points, its $what $(length(v)) values"))
        out[c] = v
    end
    out
end

# read W, the rates and ydot from a velocity's named tuple for `problem`, checked
Base.@nospecializeinfer function _read_aux(problem, @nospecialize(r), states::Vector{ThetaLState})
    nf = length(problem.field_names); bulk = problem.bulk
    nf == 0 && bulk === nothing && return _NO_AUX
    get_(k) = r isa NamedTuple && haskey(r, k) ? getproperty(r, k) : nothing
    Wr = get_(:W); W = Wr === nothing ? nothing : _checked_curve_vectors(Wr, states, "W")
    need_W = any(f -> f.transport !== :none, problem.fields)
    need_W && W === nothing &&
        throw(ArgumentError("a field is transported along the curves: the velocity must return W (`(; U, W, …)`)"))
    rr = get_(:rates); rates = Vector{Union{Nothing,Vector{Vector{Float64}}}}(nothing, nf)
    if rr !== nothing
        rr isa NamedTuple || throw(ArgumentError("`rates` must be a named tuple, one entry per field"))
        for k in keys(rr)
            i = findfirst(==(k), problem.field_names)
            i === nothing && throw(ArgumentError("the velocity returned a rate for $k, which is not a field"))
            rates[i] = _checked_curve_vectors(getproperty(rr, k), states, "rate of $k")
        end
    end
    ydot = nothing
    if bulk !== nothing
        yr = get_(:ydot)
        yr === nothing && throw(ArgumentError("the problem has a bulk state: the velocity must return ydot"))
        ydot = Vector{Float64}(yr)
        length(ydot) == length(bulk.y0) ||
            throw(DimensionMismatch("ydot has $(length(ydot)) components, the bulk state $(length(bulk.y0))"))
    end
    _Aux(W, rates, ydot)
end

# a periodic grid function resampled spectrally from N to M points (zero padding or truncation; the Nyquist mode
# split or folded as in double_grid and halve_grid)
function _resample(f::Vector{Float64}, M::Int)
    N = length(f); M == N && return copy(f)
    h = _rfft(f); nN = N ÷ 2 + 1; nM = M ÷ 2 + 1; g = zeros(ComplexF64, nM)
    for j in 1:min(nN, nM)
        g[j] = h[j]
    end
    if M > N && iseven(N)
        g[nN] = h[nN] / 2                              # the Nyquist mode of N split between ±N/2
    elseif M < N && iseven(M)
        g[nM] = 2real(h[nM])                           # ±M/2 folded into the new Nyquist mode
    end
    F = _irfft!(_work(M), Vector{Float64}(undef, M), g); s = M / N
    for j in eachindex(F)
        F[j] *= s
    end
    F
end
