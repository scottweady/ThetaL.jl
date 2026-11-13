# Periodic grid functions on N equispaced points: the FFT wavenumbers, the spectral derivative and antiderivative,
# trigonometric interpolation to 2N points, and the tail of a spectrum. Real grid functions go through half
# spectra (real FFTs); curves, as x + iy, through complex FFTs.
#
# Compile time: the package is written in loops on concrete vectors (Vector{Float64}, Vector{ComplexF64}). Each
# fused broadcast or comprehension is a type of its own, compiled on first use through the generic broadcast and
# collect machinery, and those made 90 % of the tests' time; a loop
# compiles with its function. The FFTs go through `_fft`/`_ifft`, which fix the types, and the in-place plans and
# scratch vectors of work.jl, so that a function allocates only what it returns.

_wavenumber(j::Int, N::Int) = (k = j - 1; k <= N ÷ 2 ? k : k - N)   # of the j-th FFT entry (Nyquist at N/2)
_nyquist(k::Int, N::Int) = abs(k) == N ÷ 2                           # the mode dropped by derivatives

"""FFT wavenumbers 0, 1, …, N/2, −N/2+1, …, −1."""
function _kvec(N::Int)
    k = Vector{Int}(undef, N)
    for j in 1:N
        k[j] = _wavenumber(j, N)
    end
    k
end
_absk(N::Int) = (k = _kvec(N); for j in 1:N; k[j] = abs(k[j]); end; k)

# Complex FFTs of length N into new vectors (the caller keeps them): for curves as x + iy, and full spectra
_fft(f::AbstractVector{<:Real}) = (w = _work(length(f)); _load_fft!(w, Vector{ComplexF64}(undef, w.N), f))
_fft(f::Vector{ComplexF64}) = _fwd!(_work(length(f)), copy(f))
function _ifft(f::Vector{ComplexF64})
    w = _work(length(f)); c = _bwd!(w, copy(f)); s = 1 / w.N
    for j in eachindex(c)
        c[j] *= s
    end
    c
end

# Real grid functions and their half spectra (n = N/2 + 1 coefficients, k = 0, …, N/2; work.jl)
_vec64(f::Vector{Float64}) = f
_vec64(f::AbstractVector{<:Real}) = Vector{Float64}(f)
"""The half spectrum of a real grid function, into a new vector."""
_rfft(f::AbstractVector{<:Real}) = (w = _work(length(f)); _rfft!(w, Vector{ComplexF64}(undef, w.n), _vec64(f)))
"""The real grid function of N points with half spectrum h, into a new vector (h kept)."""
function _irfft(h::Vector{ComplexF64}, N::Int)
    w = _work(N); c = copyto!(_take_h(w), h)
    out = _irfft!(w, Vector{Float64}(undef, N), c)
    _give!(w, c)
    out
end

"""Spectral derivative of a periodic grid function (Nyquist mode dropped); `_dα!` writes it into `out`."""
_dα(f::AbstractVector{<:Real}) = _dα!(Vector{Float64}(undef, length(f)), f)
function _dα!(out::Vector{Float64}, f::AbstractVector{<:Real})
    w = _work(length(f)); h = _rfft!(w, _take_h(w), _vec64(f))
    _dα_half!(h, w.N)
    _irfft!(w, out, h); _give!(w, h)
    out
end
# h ← ik h on a half spectrum (Nyquist mode dropped)
function _dα_half!(h::Vector{ComplexF64}, N::Int)
    for j in eachindex(h)
        k = j - 1
        h[j] = _nyquist(k, N) ? 0.0im : im * k * h[j]
    end
    h
end

"""Periodic antiderivative with value 0 at α = 0 (the input must have zero mean); `_intα!` writes it into `out`."""
_intα(f::AbstractVector{<:Real}) = _intα!(Vector{Float64}(undef, length(f)), f)
function _intα!(out::Vector{Float64}, f::AbstractVector{<:Real})
    w = _work(length(f)); h = _rfft!(w, _take_h(w), _vec64(f)); N = w.N
    for j in eachindex(h)
        k = j - 1
        h[j] = k == 0 || _nyquist(k, N) ? 0.0im : h[j] / (im * k)
    end
    _irfft!(w, out, h); _give!(w, h)
    F1 = out[1]
    for j in 1:N
        out[j] -= F1
    end
    out
end

"""Grid function with Fourier coefficients multiplied by the weights `wt` (a filter on the half spectrum, or on the
full one in FFT order); `_filtered!` writes it into `out`."""
_filtered(f::AbstractVector{<:Real}, wt::Vector{Float64}) = _filtered!(Vector{Float64}(undef, length(f)), f, wt)
function _filtered!(out::Vector{Float64}, f::AbstractVector{<:Real}, wt::Vector{Float64})
    w = _work(length(f)); h = _rfft!(w, _take_h(w), _vec64(f))
    for j in eachindex(h)
        h[j] *= wt[j]
    end
    _irfft!(w, out, h); _give!(w, h)
    out
end

# tail of a spectrum: the largest of |c| in the band M/2 ± M/8, relative to the largest of all (`rel = true`) or as
# it is
function _tail(c::Vector{ComplexF64}, rel::Bool)
    n = length(c); top = 0.0; band = 0.0
    for j in 1:n
        a = abs(c[j]); top = max(top, a)
        n ÷ 2 - n ÷ 8 + 1 <= j <= n ÷ 2 + n ÷ 8 + 1 && (band = max(band, a))
    end
    rel ? band / top : band
end
_rel_tail(c::Vector{ComplexF64}) = _tail(c, true)
# the same on the half spectrum of a real function of N points (the band's other half mirrors this one)
function _half_tail(h::Vector{ComplexF64}, N::Int)
    band = 0.0
    for j in N÷2-N÷8+1:length(h)
        band = max(band, abs(h[j]))
    end
    band
end

"""
    double_grid(f) / double_grid(X)

Trigonometric interpolation of a periodic grid function from N to 2N points (zero padding, the Nyquist mode split
in two), or of the φ of a state.
"""
double_grid(f::AbstractVector{<:Real}) = _double_grid(_rfft(f), length(f))
# from the half spectrum h of N points (h kept)
function _double_grid(h::Vector{ComplexF64}, N::Int)
    g = zeros(ComplexF64, N + 1)                 # the half spectrum on 2N points
    for j in 1:N÷2
        g[j] = h[j]
    end
    g[N÷2+1] = h[N÷2+1] / 2                      # the Nyquist mode split between ±N/2
    F = _irfft!(_work(2N), Vector{Float64}(undef, 2N), g)
    for j in eachindex(F)
        F[j] *= 2
    end
    F
end

"""
    halve_grid(f) / halve_grid(X)

Spectral truncation of a periodic grid function from N to N/2 points: the modes |k| < N/4 kept, the two modes
±N/4 folded into the new Nyquist mode, the rest dropped. Or of the φ of a state. The inverse of `double_grid` for a
function resolved on N/2 points.
"""
halve_grid(f::AbstractVector{<:Real}) = _halve_grid(_rfft(f), length(f))
# from the half spectrum h of N points (h kept)
function _halve_grid(h::Vector{ComplexF64}, N::Int)
    M = N ÷ 2
    M >= 2 || throw(ArgumentError("halve_grid needs at least 4 points, got $N"))
    g = Vector{ComplexF64}(undef, M ÷ 2 + 1)     # the half spectrum on M points
    for j in 1:M÷2                               # k = 0, 1, …, M/2 − 1
        g[j] = h[j]
    end
    g[M÷2+1] = 2real(h[M÷2+1])                   # k = ±M/2 folded into the new Nyquist mode
    F = _irfft!(_work(M), Vector{Float64}(undef, M), g)
    for j in eachindex(F)
        F[j] /= 2
    end
    F
end

# The trigonometric interpolant of a periodic grid function, as a function of α (the Nyquist mode split in two, so
# that it is real and passes through the samples)
struct _TrigInterpolant
    c::Vector{ComplexF64}
    k::Vector{Float64}
end
function _trig_interpolant(f::AbstractVector{<:Real})
    c = _fft(f); N = length(c)
    for j in 1:N
        c[j] /= N
    end
    _TrigInterpolant(c, Float64.(_kvec(N)))
end
# (callable structs in the short form, forwarding to a function: the long form `function (p::T)(…)` costs a
# second of parser compilation in a fresh Julia 1.12 session, paid by every precompilation)
(p::_TrigInterpolant)(a::Real) = _interpolant_at(p, Float64(a))
function _interpolant_at(p::_TrigInterpolant, a::Float64)
    v = 0.0; N = length(p.c)
    @inbounds for j in 1:N
        w = _nyquist(Int(p.k[j]), N) ? real(p.c[j]) * cos(p.k[j] * a) : real(p.c[j] * cis(p.k[j] * a))
        v += w
    end
    v
end
