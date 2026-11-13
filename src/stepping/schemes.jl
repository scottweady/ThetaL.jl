# One-step schemes. A scheme is a value carrying its settings (the spectral filter) and defining its weights; the
# generic step (src/stepping/step.jl) does everything else. Per curve and Fourier mode, with z = σh,
#
#   stage i:  φ̂⁽ⁱ⁾ = e^{cᵢ z} φ̂ + h Σⱼ Aᵢⱼ(z) Rⱼ          end:  φ̂ₙ₊₁ = e^{z} φ̂ + h Σᵢ bᵢ(z) Rᵢ
#
# with Rⱼ = θ_t − σφ̂ at stage j (the explicit remainder). L uses the same weights at z = 0, and x₀ the rigid motion
# averaged with them (exact when the motion is constant over the step). An embedded b̂ gives an error estimate.
#
# To add a scheme: a struct with one field `filter::SpectralFilter`, a keyword constructor, and methods for
# `stage_nodes`, `scheme_weights` and (optionally) `has_error_estimate`. Nothing else changes.

"""Supertype of the one-step schemes of `thetaL_step` (`ETDRK2`, `ExpEuler`)."""
abstract type OneStepScheme end

"""
    ETDRK2(; filter_order = 36, filter_cutoff = 0)

Exponential Runge–Kutta of order 2 (Hochbruck & Ostermann 2005, c₂ = 1/2), the default: two velocity evaluations
per step (the one at the end is the next step's first stage). Embedded exponential Euler: an error estimate
available after the midpoint stage, before the end evaluation. Keeps a state at rest fixed with a symbol present;
with no symbol it is the explicit midpoint rule.
"""
struct ETDRK2 <: OneStepScheme
    filter::SpectralFilter
end

"""
    ExpEuler(; filter_order = 36, filter_cutoff = 0)

Exponential Euler, order 1: one velocity evaluation per step. For checks and debugging; with no symbol it is the
explicit Euler method.
"""
struct ExpEuler <: OneStepScheme
    filter::SpectralFilter
end

# Lawson's integrating-factor midpoint rule (order 2, two evaluations per step), kept as a reference for
# tests. Not exported: it does not keep a state at rest fixed under a symbol.
struct Lawson2 <: OneStepScheme
    filter::SpectralFilter
end

for S in (:ETDRK2, :ExpEuler, :Lawson2)
    @eval $S(; filter_order::Int=36, filter_cutoff::Union{Int,Symbol}=0) =
        $S(SpectralFilter(order=filter_order, cutoff=filter_cutoff))
end

"""The same scheme with another spectral filter."""
with_filter(scheme::S, filter::SpectralFilter) where {S<:OneStepScheme} = S(filter)

"""The stage nodes cᵢ of a scheme (fractions of the step; c₁ = 0, the start of the step)."""
stage_nodes(::ETDRK2) = (0.0, 0.5)
stage_nodes(::ExpEuler) = (0.0,)
stage_nodes(::Lawson2) = (0.0, 0.5)

"""Whether the scheme has an embedded solution for an error estimate."""
has_error_estimate(::OneStepScheme) = false
has_error_estimate(::ETDRK2) = true

"""
    scheme_weights(scheme, z) -> (A, b, b̂)

The weights of one mode at z = σh (σ the symbol's rate): `A[i][j]` for stage i from stage j < i, `b[i]` for the
end, `b̂[i]` for the embedded solution (equal to `b` when there is none). Tuples, so nothing is allocated.
"""
function scheme_weights(::ETDRK2, z::Float64)
    p1h, _ = _ϕ12(z / 2); p1, p2 = _ϕ12(z)
    ((), (0.5 * p1h,)), (p1 - 2p2, 2p2), (p1, 0.0)
end
function scheme_weights(::ExpEuler, z::Float64)
    p1, _ = _ϕ12(z)
    ((),), (p1,), (p1,)
end
function scheme_weights(::Lawson2, z::Float64)
    eh = exp(z / 2)
    ((), (0.5 * eh,)), (0.0, eh), (0.0, eh)
end

"""Evaluations of the velocity per step: one per stage after the first (the first reuses the previous end), and
the end."""
evaluations_per_step(scheme::OneStepScheme) = length(stage_nodes(scheme))

# ϕ₁(z) = (eᶻ − 1)/z and ϕ₂(z) = (eᶻ − 1 − z)/z², by their series near 0 (the last term kept is z⁷)
function _ϕ12(z::T) where {T<:Union{Float64,ComplexF64}}
    if abs(z) < 0.05
        p1 = one(T); p2 = one(T) / 2; zk = one(T); f1 = 1.0; f2 = 2.0
        for k in 1:7
            zk *= z; f1 *= k + 1; f2 *= k + 2
            p1 += zk / f1; p2 += zk / f2
        end
        return p1, p2
    end
    e = expm1(z)
    e / z, (e - z) / z^2
end
