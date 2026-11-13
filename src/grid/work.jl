# Work arrays of the grid functions. For each N: the FFT plans (complex in place, for curves as x + iy; real to
# half-spectrum and back, for real grid functions), pools of scratch vectors, |k| and the filter weights on the
# half spectrum, kept per task (threads never share them) and reused across calls, so that the step and the
# geometry allocate only what they return (before, half of a step at N = 8192 went to allocating FFT temporaries).
#
# Spectra of real grid functions are half spectra: n = N/2 + 1 coefficients, k = 0, 1, …, N/2 (the Nyquist mode
# last), unnormalized as FFTW's rfft. A kernel takes a vector with `_take_c` (complex, N), `_take_h` (half
# spectrum, n), `_take_r` (real, N) or `_take_hr` (real, n), contents undefined, and hands it back with
# `_give!`; a vector it returns to its caller is never taken from a pool. A vector not handed back (an
# exception in between) is collected as usual.

const _FwdPlan = FFTW.cFFTWPlan{ComplexF64,-1,true,1,Tuple{Int}}
const _BwdPlan = FFTW.cFFTWPlan{ComplexF64,1,true,1,Tuple{Int}}
const _RFwdPlan = FFTW.rFFTWPlan{Float64,-1,false,1,Tuple{Int}}
const _RBwdPlan = FFTW.rFFTWPlan{ComplexF64,1,false,1,UnitRange{Int}}

struct _GridWork
    N::Int
    n::Int                                              # N/2 + 1, the length of a half spectrum
    fwd::_FwdPlan                                       # complex, in place, forward
    bwd::_BwdPlan                                       # complex, in place, backward, unnormalized (bfft)
    rfwd::_RFwdPlan                                     # real → half spectrum
    rbwd::_RBwdPlan                                     # half spectrum → real, unnormalized (brfft; destroys its input)
    complex::Vector{Vector{ComplexF64}}
    half::Vector{Vector{ComplexF64}}
    real::Vector{Vector{Float64}}
    halfreal::Vector{Vector{Float64}}                   # real vectors on the half spectrum (weights per mode)
    absk::Vector{Int}                                   # k = 0, 1, …, N/2: the wavenumbers of a half spectrum
    filters::Dict{Tuple{Int,Int},Vector{Float64}}       # filter weights on the half spectrum, by (order, cutoff)
    zero_rates::Vector{Float64}                         # zeros on the half spectrum (no symbol); read only
end
function _GridWork(N::Int)
    flags = FFTW.ESTIMATE | FFTW.UNALIGNED              # the plans apply to any vectors of the right lengths
    x = zeros(ComplexF64, N); r = zeros(N); h = zeros(ComplexF64, N ÷ 2 + 1)
    _GridWork(N, N ÷ 2 + 1, FFTW.plan_fft!(x; flags=flags), FFTW.plan_bfft!(x; flags=flags),
              FFTW.plan_rfft(r; flags=flags), FFTW.plan_brfft(h, N; flags=flags), Vector{ComplexF64}[],
              Vector{ComplexF64}[], Vector{Float64}[], Vector{Float64}[], collect(0:N÷2),
              Dict{Tuple{Int,Int},Vector{Float64}}(), zeros(N ÷ 2 + 1))
end

# the work arrays of this task for N points
function _work(N::Int)
    d = get!(Dict{Int,_GridWork}, task_local_storage(), :ThetaL_grid_work)::Dict{Int,_GridWork}
    w = get(d, N, nothing)
    w === nothing ? (d[N] = _GridWork(N)) : w
end

_take_c(w::_GridWork) = isempty(w.complex) ? Vector{ComplexF64}(undef, w.N) : pop!(w.complex)
_take_h(w::_GridWork) = isempty(w.half) ? Vector{ComplexF64}(undef, w.n) : pop!(w.half)
_take_r(w::_GridWork) = isempty(w.real) ? Vector{Float64}(undef, w.N) : pop!(w.real)
_take_hr(w::_GridWork) = isempty(w.halfreal) ? Vector{Float64}(undef, w.n) : pop!(w.halfreal)
_give!(w::_GridWork, v::Vector{ComplexF64}) = (push!(length(v) == w.N ? w.complex : w.half, v); nothing)
_give!(w::_GridWork, v::Vector{Float64}) = (push!(length(v) == w.N ? w.real : w.halfreal, v); nothing)
_give!(w::_GridWork, v, vs...) = (_give!(w, v); _give!(w, vs...))

# in place: the forward FFT, and the backward one without the 1/N
_fwd!(w::_GridWork, c::Vector{ComplexF64}) = (w.fwd * c; c)
_bwd!(w::_GridWork, c::Vector{ComplexF64}) = (w.bwd * c; c)

# f (real) into c as complex numbers, transformed (a full spectrum, for `_fft` of real data)
function _load_fft!(w::_GridWork, c::Vector{ComplexF64}, f::AbstractVector{<:Real})
    for j in 1:w.N
        c[j] = f[j]
    end
    _fwd!(w, c)
end

# real grid function ↔ half spectrum: f into h (`_rfft!`), and the inverse of h into `out` with the 1/N (`_irfft!`;
# h is overwritten)
_rfft!(w::_GridWork, h::Vector{ComplexF64}, f::Vector{Float64}) = mul!(h, w.rfwd, f)
function _irfft!(w::_GridWork, out::Vector{Float64}, h::Vector{ComplexF64})
    mul!(out, w.rbwd, h); s = 1 / w.N
    for j in 1:w.N
        out[j] *= s
    end
    out
end
