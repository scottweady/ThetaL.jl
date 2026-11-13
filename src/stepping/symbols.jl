# The stiff symbol: the growth rate σ(q) of a small-scale mode of a curve, θ_t ≈ σ φ̂ at wavenumber q = 2π|k|/L
# along the curve (σ < 0 damps), which the exponential step treats exactly. Given by the caller as a function of q,
# `σ(q)`, or of q and the curve, `σ(q, X)`; frozen per curve at the start of each step.

# the growth rates σ_k of one curve X at the wavenumbers k (`nothing`: no symbol, zero rates), checked; the caller
# only reads them
Base.@nospecializeinfer function _curve_rates(symbol, X::ThetaLState, k::Vector{Int})
    @nospecialize symbol
    if symbol === nothing                                   # zeros; on the half spectrum shared (read only)
        return length(k) == length(X.φ) ÷ 2 + 1 ? _work(length(X.φ)).zero_rates : zeros(length(k))
    end
    withX = applicable(symbol, 1.0, X)
    withX || applicable(symbol, 1.0) ||
        throw(ArgumentError("the symbol must be a function of the wavenumber, σ(q), or of it and the curve, σ(q, X)"))
    a = 2π / X.L; σ = Vector{Float64}(undef, length(k)); positive = false
    for j in eachindex(k)
        q = a * abs(k[j])
        r = withX ? symbol(q, X) : symbol(q)
        r isa Real && isfinite(r) ||
            throw(ArgumentError("the symbol returned $r at q = $q: a finite real number is needed"))
        σ[j] = r; positive |= r > 0
    end
    positive && @warn "the symbol has positive rates σ > 0: the exponential step amplifies those modes" maxlog = 1
    σ
end
