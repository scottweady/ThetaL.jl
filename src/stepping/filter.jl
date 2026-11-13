# The spectral filter applied in every step (`SpectralFilter`), and a cutoff above the noise that a computed velocity
# carries (`noise_cutoff`).

"""
    SpectralFilter(; order = 36, cutoff = 0)

The exponential filter of the θ–L step: mode k is multiplied by exp(−36 (|k|/kc)^order), and modes |k| ≥ kc are
removed. It is applied to φ after every stage and to the deformation part of U before it is differentiated.

- `order`: the steepness; `0` turns the filter off.
- `cutoff`: kc as an integer; `0` means kc = N/2 (the Nyquist mode); `:auto` means "set it from `noise_cutoff`
  after the first velocity evaluation" (the integrator does this once and keeps the integer).
"""
Base.@kwdef struct SpectralFilter
    order::Int = 36
    cutoff::Union{Int,Symbol} = 0
end

# the filter weights on the FFT grid of N points (kc > 0: the cutoff, otherwise N/2)
function _filter(N::Int, order::Int, kc::Int=0)
    f = ones(N)
    order == 0 && return f
    km = kc > 0 ? min(kc, N ÷ 2) : N ÷ 2
    for j in 1:N
        k = abs(_wavenumber(j, N))
        f[j] = k >= km ? 0.0 : exp(-36 * (k / km)^order)
    end
    f
end
function _filter(N::Int, filter::SpectralFilter)
    filter.cutoff isa Int || throw(ArgumentError("the filter cutoff :auto has not been resolved yet"))
    _filter(N, filter.order, filter.cutoff)
end
# the weights on the half spectrum (k = 0, …, N/2), computed once per task and N (work.jl); read only
function _cached_filter(N::Int, filter::SpectralFilter)
    filter.cutoff isa Int || throw(ArgumentError("the filter cutoff :auto has not been resolved yet"))
    get!(() -> _filter(N, filter.order, filter.cutoff)[1:N÷2+1], _work(N).filters, (filter.order, filter.cutoff))
end

"""
    noise_cutoff(states, U, symbol) -> kc

A filter cutoff above the noise that the normal velocities `U` of the `states` carry. An error δÛ_k in U about flat
in k, and all that U holds at high k for a smooth start, is measured as ε = max |Û_k|/N over k ≥ N/4. It enters φ
through θ_t ∋ (2π/L)(ik)Û_k, against which the step damps at the rate −σ_k from `symbol`, so the noise settles at
φ̂_k ≈ (2π/L) k ε / |σ_k|. The cutoff is placed so the filter keeps 99 % of the mode where the φ spectrum meets this
floor (σ = 0.99 at 0.8 kc), and is at least 32, so curves that start as circles can deform.
"""
Base.@nospecializeinfer function noise_cutoff(states::AbstractVector{ThetaLState}, U::AbstractVector, symbol)
    @nospecialize symbol
    kx = 1
    for (X, Uc) in zip(states, U)
        N = length(X.φ); n = N ÷ 2 + 1
        cU = _rfft(Uc); cφ = _phihat(X)                    # half spectra: k = 0, …, N/2
        ε = 0.0
        for j in N÷4+1:n
            ε = max(ε, abs(cU[j]) / N)
        end
        r = -_curve_rates(symbol, X, collect(0:N÷2))                  # the damping rates −σ_k, k = 0, …, N/2
        m = n                                                 # φ̂_k is below the floor for every k ≥ m
        while m > 2 && abs(cφ[m]) / N < (2π / X.L) * (m - 1) * ε / r[m]   # (entry m is k = m − 1)
            m -= 1
        end
        kx = max(kx, m)
    end
    max(ceil(Int, kx / 0.8), 32)
end
