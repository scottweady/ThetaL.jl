# The stiffness of a velocity measured against its symbol (src/problem/symbol_check.jl), no step.
#   measure_symbol: on the unit circle the measured rates are the exact linearization of the step's θ_t (curve
#   shortening 2 − q², surface diffusion q² − q⁴, a Hilbert transform −c(q − 1/q)); on a non-round curve close to
#   the leading term; local windows see a stiffness that varies along the curve, which the global mode averages.
#   check_symbol: the stability bound per mode, the right symbol passes, a weak one or none is flagged.
#   julia -t 8 --project -e 'using Pkg; Pkg.test(test_args=["symbol_check"])'   (from the package root)
using ThetaL, LinearAlgebra, Test
using ThetaL: Vec2, _stable_step, Lawson2
import ThetaL.FFTW
include(joinpath(@__DIR__, "common", "curves.jl"))

circle(N) = ThetaLState(zeros(N), 2π, Vec2(0, -1))                  # the unit circle: q = k
function hilbert(f)                                                  # symbol −i sgn k
    N = length(f); c = FFTW.fft(f); k = FFTW.fftfreq(N, N)
    real(FFTW.ifft(-im .* sign.(k) .* c))
end
const MCF = (Xs, ev; warm) -> [-curvature(X) for X in Xs]
const SD = (Xs, ev; warm) -> [arclength_derivative(X, arclength_derivative(X, curvature(X))) for X in Xs]
hilbert_velocity(c) = (Xs, ev; warm) -> [c .* hilbert(X.φ) for X in Xs]

@testset "measured rates on a circle: the exact linearization" begin
    c = 0.7
    for (U, σ, exact) in ((MCF, q -> -q^2, q -> 2 - q^2), (SD, q -> -q^4, q -> q^2 - q^4),
                          (hilbert_velocity(c), q -> -c * q, q -> -c * (q - 1 / q)))
        r = measure_symbol(CurveProblem([circle(128)], U; symbol=σ))
        @test r.k == [2, 3, 6, 11, 18, 32] && r.q ≈ r.k
        @test all(isapprox.(r.σ, exact.(r.q); rtol=1e-6))
        @test maximum(abs, imag.(r.λ)) < 1e-6 * maximum(abs, r.σ)
        @test r.given == σ.(r.q)
        @test abs(r.σ[end] / r.given[end] - 1) < 3e-3                 # the leading term at q = 32
    end
end

@testset "measured rates on a non-round curve" begin
    X = smooth_state(128); c = 0.7
    for (U, σ) in ((MCF, q -> -q^2), (SD, q -> -q^4), (hilbert_velocity(c), q -> -c * q))
        r = measure_symbol(CurveProblem([X], U; symbol=σ); k=[8, 16, 32])
        @test all(abs.(r.σ ./ r.given .- 1) .< [0.05, 0.015, 0.006])
    end
    # the other curve held; the probes are told :probe and the context is passed along
    kinds = Symbol[]; warms = Any[]
    U2 = (Xs, ev; warm) -> (push!(kinds, ev.kind); push!(warms, warm);
                            (; U=[-curvature(Xs[1]), zeros(length(Xs[2]))], ctx=length(warms)))
    r = measure_symbol(CurveProblem([circle(64), X], U2; symbol=q -> -q^2); curve=1, k=[4], warm=:w)
    @test kinds == fill(:probe, 4) && warms == [:w, 1, 2, 3]
    @test r.σ[1] ≈ 2 - 16 rtol = 1e-6
    r = measure_symbol(CurveProblem([circle(64), X], U2; symbol=q -> -q^2); curve=2, k=[4])
    @test abs(r.σ[1]) < 1e-12                                        # curve 2 does not move
end

@testset "local windows: a stiffness that varies along the curve" begin
    # U = −a(α)κ, a = 1 + 0.8 cos α on the circle: the global mode sees the mean a = 1, the windows a(α₀)
    a = [1 + 0.8cos(2π * (j - 1) / 128) for j in 1:128]
    U = (Xs, ev; warm) -> [-a .* curvature(Xs[1])]
    p = CurveProblem([circle(128)], U; symbol=q -> -q^2)
    r1 = measure_symbol(p; k=[32]); r4 = measure_symbol(p; k=[32], windows=4)
    @test r1.σ[1] ≈ 2 - 32^2 rtol = 1e-3
    @test size(r4.λ) == (1, 4)
    @test real(r4.λ[1, 1]) / -32^2 > 1.7 && real(r4.λ[1, 3]) / -32^2 < 0.3      # α₀ = 0 and π
    @test r4.σ[1] == minimum(real, r4.λ)
    @test check_symbol(p; k=[32]).ok
    @test !(@test_logs (:warn,) check_symbol(p; k=[32], windows=4)).ok
    @test check_symbol(CurveProblem([circle(128)], U; symbol=q -> -1.8q^2); k=[32], windows=4).ok
end

@testset "the stability bound of a mode" begin
    # σ_g treated exactly, λ − σ_g explicitly: at large steps ETDRK2 needs 0 ≤ λ/σ_g ≤ 3/2, ExpEuler 0 ≤ λ/σ_g ≤ 2
    σg = -100.0
    for (ρ, rk2, ee) in ((1.0, true, true), (0.0, true, true), (1.45, true, true), (1.55, false, true),
                         (1.95, false, true), (2.05, false, false), (-0.5, true, true))
        @test isinf(_stable_step(ETDRK2(), σg, complex(ρ * σg))) == rk2
        @test isinf(_stable_step(ExpEuler(), σg, complex(ρ * σg))) == ee
    end
    # no symbol: ETDRK2 is the explicit midpoint rule, stable for h|λ| ≤ 2
    @test _stable_step(ETDRK2(), 0.0, complex(-100.0)) ≈ 0.02 rtol = 1e-6
    # a rate that truly grows may grow as fast as it does
    @test isinf(_stable_step(ETDRK2(), 0.0, complex(3.0)))
    @test isinf(_stable_step(ETDRK2(), 0.0, 0.0im))
end

@testset "check_symbol" begin
    X = circle(128)
    @test check_symbol(CurveProblem([X], MCF; symbol=q -> -q^2)).ok
    @test check_symbol(CurveProblem([smooth_state(128)], SD; symbol=q -> -q^4)).ok
    @test check_symbol(CurveProblem([X], hilbert_velocity(0.7); symbol=q -> -0.7q)).ok
    for σ in (q -> -q^2 / 10, nothing)
        p = CurveProblem([X], MCF; symbol=σ)
        r = @test_logs (:warn, r"does not bound") check_symbol(p)
        @test !r.ok && r.dt_limit < 0.003 && r.dt_limit == minimum(r.curves[1].dt_max)
        @test all(isfinite, r.curves[1].dt_max)
        @test issorted(r.curves[1].dt_max; rev=true)                  # the limit set by the highest q
        # a step below the limit passes, one above fails
        @test check_symbol(p; dt=r.dt_limit / 2, warn=false).ok
        @test !check_symbol(p; dt=2r.dt_limit, warn=false).ok
    end
    # λ/σ_g = 1.8: beyond ETDRK2's 3/2, within ExpEuler's 2
    p = CurveProblem([X], MCF; symbol=q -> -q^2 / 1.8)
    @test !check_symbol(p; warn=false).ok && check_symbol(p; scheme=ExpEuler()).ok
    # the step confirms the limit to a few percent: U = 2π/L − κ (a circle at rest) under −q²/10, no filter, every
    # mode up to N/2 − 1 measured; seeded at small scales, which decay below the limit and grow above it
    N = 64; Y = ThetaLState([0.2sin(2α) + 1e-4 * sum(cos(k * α) for k in 20:N÷2-1) for α in 2π .* (0:N-1) ./ N],
                            20π, Vec2(0.0, -10.0))
    σw = q -> -q^2 / 10; V = (Z, t) -> 2π / Z.L .- curvature(Z)
    lim = check_symbol(CurveProblem([Y], (Zs, ev; warm) -> [V(Z, 0.0) for Z in Zs]; symbol=σw); k=2:N÷2-1,
                       warn=false).dt_limit
    tail(Z) = maximum(abs, ThetaL._rfft(Z.φ)[21:N÷2]) / N
    steps(dt) = try
        tail(evolve(Y, V, (0, 200dt); dt, symbol=σw, scheme=ETDRK2(filter_order=0), Nmax=N).curves[end])
    catch
        Inf                                                          # blown up: a step rejected
    end
    @test steps(0.95lim) < 1e-6 * tail(Y)
    @test steps(1.05lim) > 100 * tail(Y)
    # every curve, and the integrator's method (its scheme, time and states)
    p = CurveProblem([X, smooth_state(64)], (Xs, ev; warm) -> [-curvature(Z) for Z in Xs]; symbol=q -> -q^2)
    r = check_symbol(p; k=[4, 8])
    @test r.ok && length(r.curves) == 2 && r.curves[2].curve == 2
    integ = init(p, ExpEuler())
    @test check_symbol(integ; curves=[2], k=[8]).ok
    @test measure_symbol(integ; curve=2, k=[8]).σ ≈ r.curves[2].σ[2:2]
    @test_throws ArgumentError measure_symbol(p; curve=3)
    @test_throws ArgumentError measure_symbol(p; k=[64])
end
