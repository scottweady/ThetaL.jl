# The step (src/stepping/), through `evolve` and `solve`, under prescribed velocities and symbols; no solve.
#   Kinematics: a uniform normal velocity on a circle (exact); the rigid part of a velocity; the frame velocity.
#   Schemes (ETDRK2, ExpEuler; Lawson2 as a reference): rigid motions exact at any Δt with and without a symbol;
#   no symbol gives the classical method; order under a linear flow; order from a relaxed state (surface
#   diffusion); a manufactured deforming solution under wrong symbols; the embedded estimate; a pure stiff term;
#   the rigid split on a perturbed shrinking near-circle. The filter and its noise cutoff (no step).
#   julia -t 8 --project -e 'using Pkg; Pkg.test(test_args=["stepping"])'   (from the package root)
using ThetaL, LinearAlgebra, Printf, Test
using ThetaL: Vec2, _kvec, _filter, Lawson2
import ThetaL.FFTW
include(joinpath(@__DIR__, "common", "curves.jl"))

const SCHEMES = (ETDRK2(), ExpEuler())

@testset "uniform normal velocity on a circle" begin
    N, c = 64, 0.7
    X0 = ThetaLState(zeros(N), 2π, Vec2(0, -1))                      # unit circle; x(0) = (0, −1), τ(0) = (1, 0)
    for scheme in SCHEMES
        sol = evolve(X0, (X, t) -> fill(c, N), (0, 0.5); dt=0.1, symbol=q -> -q, scheme)
        X = sol.curves[end]; R = 1 + 0.5c
        @test sol.steps == 5
        @test X.L ≈ 2π * R rtol = 1e-14
        @test maximum(abs, X.φ) < 1e-14
        @test X.x0 ≈ Vec2(0, -R) atol = 1e-14
    end
end

@testset "the rigid part of a velocity, and the frame velocity" begin
    # a rigid motion plus a deformation orthogonal to it: the fit returns the motion, and the rest
    X = smooth_state(128); x = curve_points(X); θ = tangent_angles(X); xc = sum(x) / 128
    m0 = RigidMotion(Vec2(0.3, -0.2), 0.7, xc)
    Ur = [dot(rigid_velocity(m0, p), Vec2(sin(q), -cos(q))) for (p, q) in zip(x, θ)]
    m, Ud = rigid_part(X, Ur)
    @test norm(m.translation - m0.translation) < 1e-13 && abs(m.rotation_rate - 0.7) < 1e-13 && maximum(abs, Ud) < 1e-13
    m, Ud = rigid_part(X, Ur; mode=:translation)
    @test m.rotation_rate == 0 && maximum(abs, Ud) > 0.1
    m, Ud = rigid_part(X, Ur; mode=:none)
    @test m == ThetaL.NO_RIGID_MOTION && Ud == Ur
    # a circle: the rotation gives no normal velocity and is left out
    m, Ud = rigid_part(ThetaLState(zeros(64), 2π, Vec2(0, -1)), fill(0.5, 64))
    @test m.rotation_rate == 0 && norm(m.translation) < 1e-15 && maximum(abs.(Ud .- 0.5)) < 1e-15
    # a rigid motion: the tracked points move with it
    @test maximum(norm.(frame_velocity(X, Ur) .- rigid_velocity.(Ref(m0), x))) < 1e-12
end

@testset "rigid motions are exact at any Δt" begin
    # a smooth curve turning about a point outside it, ten steps a turn, with and without a symbol; a translation
    N, ω = 128, 2π
    X0 = smooth_state(N); z0 = curve_function(X0)
    for scheme in (SCHEMES..., Lawson2()), symbol in (nothing, q -> -q^3)
        scheme isa Lawson2 && symbol !== nothing && continue          # Lawson does not keep a state at rest
        sol = evolve(X0, field(x -> ω * Vec2(-x[2] + 0.5, x[1] - 2)), (0, 1); dt=0.1, symbol, scheme)
        e = curve_distance(z0, curve_points(sol.curves[end]))
        @printf("rigid rotation, %s, symbol %s: distance %.1e\n", nameof(typeof(scheme)), symbol === nothing ? "none" : "k³", e)
        @test sol.steps == 10 && e < 1e-12 && abs(sol.curves[end].L / X0.L - 1) < 1e-13
    end
    sol = evolve(X0, field(x -> Vec2(1.5, -0.5)), (0, 2); dt=0.5)
    @test norm(sol.curves[end].x0 - (X0.x0 + Vec2(3, -1))) < 1e-13 && maximum(abs.(sol.curves[end].φ .- X0.φ)) < 1e-13
    # with the rigid part not taken out (`rigid = :none`) the rotation is an explicit transport: not exact
    turned(p) = Vec2(cos(ω * 0.5) * p[1] + sin(ω * 0.5) * p[2], -sin(ω * 0.5) * p[1] + cos(ω * 0.5) * p[2])
    sol = evolve(X0, field(x -> ω * Vec2(-x[2], x[1])), (0, 0.5); dt=0.001, rigid=:none)   # (stable only for small Δt)
    @test 1e-10 < curve_distance(z0, turned.(curve_points(sol.curves[end]))) < 1e-3
end

@testset "second order under a linear flow; the classical method with no symbol" begin
    A = [0.3 0.5; -0.4 -0.1]; T = 0.4
    X0 = smooth_state(128); z0 = curve_function(X0)
    back = x -> Vec2(exp(-T * A) * x)
    for (scheme, order) in ((ETDRK2(), 2), (ExpEuler(), 1))
        e = map((0.1, 0.05, 0.025)) do dt
            X = evolve(X0, field(x -> Vec2(A * x)), (0, T); dt, symbol=q -> -0.5q, scheme).curves[end]
            curve_distance(z0, back.(curve_points(X)))
        end
        @printf("linear flow, %s: errors %.2e %.2e %.2e, ratios %.2f %.2f\n", nameof(typeof(scheme)), e..., e[1] / e[2], e[2] / e[3])
        @test 0.85 * 2^order < e[1] / e[2] < 1.15 * 2^order && 0.85 * 2^order < e[2] / e[3] < 1.15 * 2^order
    end
    # no symbol: ETDRK2 is the explicit midpoint rule, ExpEuler explicit Euler; with U = t on a circle,
    # L_t = 2πt, the midpoint rule integrates it exactly and Euler lags by 2π Σ h·t_n
    for (scheme, L) in ((ETDRK2(), 2π * (1 + 0.3^2 / 2)), (ExpEuler(), 2π * (1 + 0.1 * (0.0 + 0.1 + 0.2))))
        X = evolve(ThetaLState(zeros(64), 2π, Vec2(0, -1)), (X, t) -> fill(t, 64), (0, 0.3); dt=0.1, scheme).curves[end]
        @test X.L ≈ L rtol = 1e-14
    end
end

@testset "surface diffusion from a relaxed state: order" begin
    # U = κ_ss with its symbol σ = −q⁴, a three-lobed curve on 64 points, first relaxed for t = 0.005 with a small
    # step (from the curve itself the first steps carry an initial layer, first order in Δt: the fast decay of its
    # small scales drives L_t, which the explicit L update follows only on Δt ≪ 1/λ); then Δt = 4, 2, 1 × 10⁻³ to
    # t = 0.05 against Δt = 2.5 × 10⁻⁵
    κss(X, t) = arclength_derivative(X, arclength_derivative(X, curvature(X)))
    sd = q -> -q^4
    X1 = evolve(from_curve(θ -> (1 + 0.1cos(3θ)) .* (cos(θ), sin(θ)), 64), κss, (0, 0.005); dt=1e-5, symbol=sd,
                Nmax=64).curves[end]
    for (scheme, order) in ((ETDRK2(), 2), (ExpEuler(), 1))
        ref = evolve(X1, κss, (0, 0.05); dt=2.5e-5, symbol=sd, Nmax=64, scheme).curves[end]
        e = map((4e-3, 2e-3, 1e-3)) do dt
            abs(evolve(X1, κss, (0, 0.05); dt, symbol=sd, Nmax=64, scheme).curves[end].L - ref.L)
        end
        @printf("surface diffusion, %s: L errors %.2e %.2e %.2e, ratios %.2f %.2f\n", nameof(typeof(scheme)), e...,
                e[1] / e[2], e[2] / e[3])
        @test 0.85 * 2^order < e[1] / e[2] < 1.15 * 2^order && 0.85 * 2^order < e[2] / e[3] < 1.15 * 2^order
    end
end

@testset "a manufactured deforming solution under wrong symbols" begin
    # an ellipse with axes a(t), b(t) turning at ω: the curve F(x, t) = 0, F = (x′/a)² + (y′/b)² − 1, x′ = R(−ωt)x,
    # moves with the normal velocity U = −F_t/|∇F| at any point. Nothing in U is stiff; the symbols given are wrong
    # on purpose (0, ×0.5 and ×2 of a first-order one), and the step must still converge at second order. (Δt is
    # kept below the limit of the explicit tangential transport of the deformation, about Δt ≲ s_α/(k_max |T|):
    # at Δt = 0.02 the runs without a damping symbol are unstable.)
    ω = 1.0; a(t) = 1.5 + 0.3sin(t); b(t) = 0.8 - 0.2sin(2t); ȧ(t) = 0.3cos(t); ḃ(t) = -0.4cos(2t)
    rot(x, t) = (cs = cos(ω * t); sn = sin(ω * t); Vec2(cs * x[1] + sn * x[2], -sn * x[1] + cs * x[2]))   # x′ = R(−ωt)x
    function level(x, t)
        p = rot(x, t); F = (p[1] / a(t))^2 + (p[2] / b(t))^2 - 1
        g = Vec2(2p[1] / a(t)^2, 2p[2] / b(t)^2)                     # ∇′F (|∇F| = |∇′F|)
        Ft = -2p[1]^2 * ȧ(t) / a(t)^3 - 2p[2]^2 * ḃ(t) / b(t)^3 + dot(g, Vec2(ω * p[2], -ω * p[1]))
        F, Ft, norm(g)
    end
    U(X, t) = [(r = level(x, t); -r[2] / r[3]) for x in curve_points(X)]
    X0 = from_curve(t -> Vec2(a(0) * cos(t), b(0) * sin(t)), 256); T = 1.0
    for scale in (0.0, 0.5, 2.0)
        e = map((0.01, 0.005, 0.0025)) do dt
            X = evolve(X0, U, (0, T); dt, symbol=q -> -scale * q, Nmax=256).curves[end]
            maximum(x -> (r = level(x, T); abs(r[1]) / r[3]), curve_points(X))  # ≈ distance to the exact curve
        end
        @printf("manufactured ellipse, symbol ×%.1f: errors %.2e %.2e %.2e, ratios %.2f %.2f\n", scale, e..., e[1] / e[2], e[2] / e[3])
        @test 3.4 < e[1] / e[2] < 4.6 && 3.4 < e[2] / e[3] < 4.6
    end
end

@testset "the embedded error estimate" begin
    # ETDRK2 against its embedded exponential Euler: O(Δt²) per step
    A = [0.3 0.5; -0.4 -0.1]
    prob = CurveProblem([smooth_state(128)], field_velocity(x -> Vec2(A * x)); symbol=q -> -0.5q)
    e = map((0.02, 0.01)) do dt
        integ = init(prob, ETDRK2()); advance!(integ, dt); integ.stats.error_estimate
    end
    @test 3.6 < e[1] / e[2] < 4.4
    integ = init(prob, ExpEuler()); advance!(integ, 0.01)
    @test isnan(integ.stats.error_estimate)
end

@testset "exponential step against a pure stiff term" begin
    # for φ = O(ε) the normal velocity Û_k = iσ_k k s_α φ̂_k/(k² − 1) gives θ_t = σ_k φ̂_k + O(ε²), which the
    # exponential step treats exactly, also at |σ|Δt > 1 where an explicit step is unstable (σ = −50q: q = |k| on
    # L = 2π)
    N, ε, dt, n, m = 64, 1e-6, 0.01, 4, 3
    symbol = q -> -50q
    function U(X, t)
        k = _kvec(length(X.φ)); sα = X.L / 2π; σ = -50.0 .* abs.(k)
        Uh = [abs(k[i]) < 2 ? 0.0im : im * σ[i] * k[i] * sα / (k[i]^2 - 1) for i in eachindex(k)] .* FFTW.fft(X.φ)
        real.(FFTW.ifft(Uh))
    end
    X0 = ThetaLState([ε * cos(m * α) for α in 2π .* (0:N-1) ./ N], 2π, Vec2(0, -1))
    sol = evolve(X0, U, (0, n * dt); dt, symbol)
    g = exp(-50.0 * m * n * dt)
    φ = sol.curves[end].φ .- sum(sol.curves[end].φ) / N
    err = maximum(abs.(φ .- g .* X0.φ)) / (ε * g)
    @printf("stiff term, |σ|Δt = %.1f: decay %.3e, relative error %.1e\n", 50.0 * m * dt, g, err)
    @test sol.steps == n && err < 1e-3
end

@testset "the rigid split on a perturbed shrinking near-circle" begin
    # curve shortening of a circle with a 10⁻⁴ perturbation: the rotation fit sees a small b and may return a large Ω,
    # but the split is exact, so the curve agrees with the one stepped without a rigid part
    X0 = from_curve(t -> (1 + 1e-4cos(3t)) .* Vec2(cos(t), sin(t)), 64)
    run(rigid) = evolve(X0, (X, t) -> -curvature(X), (0, 0.2); dt=1e-3, symbol=q -> -q^2, rigid).curves[end]
    a, b = run(:full), run(:none)
    d = curve_distance(curve_function(b), curve_points(a))
    @printf("near-circle, rigid :full against :none: distance %.1e\n", d)
    @test d < 1e-9 && abs(a.L - b.L) < 1e-12
end

@testset "a near-circle with noise in U: the spurious rotation does not move the curve" begin
    # within 1e-6 of a circle U hardly shows a rotation, so 1e-6 of noise in U fits an Ω of order one, different at
    # each stage. That turns the curve's points along it, but the curve must stay the one stepped without a rigid
    # part: x₀'s drift turns with the curve (each stage's drift in the frame of the start of the step)
    N, η = 128, 1e-6
    X0 = from_curve(ξ -> (cos(ξ) * (1 + 1e-6cos(2ξ)), sin(ξ)), N; tol=1e-13); α = 2π .* (0:N-1) ./ N
    function run(rigid)
        k = 0; Ωmax = 0.0
        noise() = (k += 1; η .* [sum(sin(m * a + k + m) / m for m in 1:4) for a in α])   # the same at every call k
        vel(Xs, ev; warm) = [-curvature(X) .+ noise() for X in Xs]
        integ = init(CurveProblem([X0], vel; symbol=q -> -q^2, rigid), ETDRK2())
        for j in 1:20
            advance!(integ, j * 0.005); Ωmax = max(Ωmax, abs(integ.rigid[1].rotation_rate))
        end
        integ.states[1], Ωmax
    end
    (a, Ω), (b, _) = run(:full), run(:none)
    d = curve_distance(curve_function(b), curve_points(a)); dc = norm(centroid(a) - centroid(b))
    @printf("near-circle with noise, :full against :none: max |Ω| %.2f, distance %.1e, centroids %.1e\n", Ω, d, dc)
    @test Ω > 0.1                                   # the case is the one meant: a large spurious rotation
    @test d < 1e-9 && dc < 1e-9 && abs(a.L - b.L) < 1e-12
end

@testset "the rigid fit from the whole velocity (U and W)" begin
    # u = V + Ω ẑ × (x − p): from U alone a circle's rotation is invisible, from (U, W) it is exact for any curve
    u(x) = Vec2(0.3, -0.2) + 0.7 * Vec2(-(x[2] - 0.4), x[1] + 1.1)
    for X in (ThetaLState(zeros(64), 2π, Vec2(0, -1)), smooth_state(128))
        U, W = velocity_components(X, u.(curve_points(X)))
        m, Ud = rigid_part(X, U; W)
        @test m.rotation_rate ≈ 0.7 atol = 1e-13
        @test maximum(norm(rigid_velocity(m, x) - u(x)) for x in curve_points(X)) < 1e-13
        @test maximum(abs, Ud) < 1e-13
        @test maximum(norm.(frame_velocity(X, U; W) .- u.(curve_points(X)))) < 1e-13
    end
end

@testset "filter" begin
    N, kc = 64, 16
    f = _filter(N, 36, kc); k = _kvec(N)
    @test f[1] == 1 && all(f[abs.(k) .>= kc] .== 0)
    @test f ≈ [abs(q) >= kc ? 0.0 : exp(-36 * (abs(q) / kc)^36) for q in k]
    @test _filter(N, 0) == ones(N) && _filter(N, 36) == _filter(N, 36, N ÷ 2)
    # one step at rest: the modes at and above the cutoff are removed, the low ones kept
    X0 = ThetaLState([1e-3 * (cos(2α) + cos(20α)) for α in 2π .* (0:N-1) ./ N], 2π, Vec2(0, -1))
    X = evolve(X0, (X, t) -> zeros(N), (0, 0.1); dt=0.1, scheme=ETDRK2(filter_cutoff=kc)).curves[end]
    @test maximum(abs.(X.φ .- [1e-3 * f[3] * cos(2α) for α in 2π .* (0:N-1) ./ N])) < 1e-16
end

@testset "noise cutoff of the filter" begin
    # |φ̂_k|/N = 10^(−k/8) against a flat noise ε = 10^(−6.05) in U, σ = −q = −k (L = 2π): the floor is ε, which the
    # spectrum meets at k = 48, so kc = ⌈49/0.8⌉ = 62; a noise above the whole spectrum gives the minimum 32
    N = 256; α = 2π .* (0:N-1) ./ N
    φ = [sum(2 * 10^(-k / 8) * cos(k * a) for k in 1:N÷2-1) for a in α]
    noise(ε) = [j == 1 ? ε * N : 0.0 for j in 1:N]
    rate = q -> -q
    @test noise_cutoff([ThetaLState(φ, 2π, Vec2(0, -1))], [noise(10^-6.05)], rate) == 62
    @test noise_cutoff([ThetaLState(φ, 2π, Vec2(0, -1))], [noise(1.0)], rate) == 32
end
