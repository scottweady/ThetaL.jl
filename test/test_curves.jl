# States of curves (src/curves/, src/grid/fourier.jl): from arclength samples and from a parametrized curve
#   (`from_curve`: an ellipse against its exact length and area, a known state through another parametrization, the
#   automatic N); the geometry of a state (moments, area, centroid, perimeter, curvature, arclength derivative,
#   velocity components); its resolution (spectral tail, doubling and halving N); the measures for callers
#   (distances, integrals, the parallel curve, area rates)
#   julia -t 8 --project -e 'using Pkg; Pkg.test(test_args=["curves"])'   (from the package root)
using ThetaL, LinearAlgebra, Printf, Test
using ThetaL: Vec2
using SpecialFunctions: ellipe
include(joinpath(@__DIR__, "common", "curves.jl"))

@testset "states from arclength samples" begin
    X = smooth_state(128)
    r = from_arclength_samples(curve_points(X))
    @test maximum(abs.(r.φ .- X.φ)) < 1e-12 && abs(r.L - X.L) < 1e-12 && r.x0 ≈ X.x0
    X = from_arclength_samples([Vec2(2cos(α), 2sin(α)) for α in 2π .* (0:31) ./ 32])   # starts at (2, 0), going up
    @test X.L ≈ 4π && maximum(abs.(X.φ .- π / 2)) < 1e-12
end

@testset "states from a parametrized curve" begin
    # the ellipse 1.5 × 0.75: on the curve, its length (complete elliptic integral) and area; the points at equal
    # arclength, so that reading them as such gives the same state
    γ = θ -> Vec2(1.5cos(θ), 0.75sin(θ))
    X = from_curve(γ, 128)
    # 6.0e-11 on the curve: the accuracy of the curve rebuilt from 128 tangent angles
    @test maximum(abs(sqrt((x[1] / 1.5)^2 + (x[2] / 0.75)^2) - 1) for x in curve_points(X)) < 1e-10
    @test X.L ≈ 6ellipe(0.75) rtol = 1e-13
    @test curve_moments(X)[1] ≈ π * 1.5 * 0.75 rtol = 1e-12
    @test X.x0 == Vec2(1.5, 0)
    r = from_arclength_samples(curve_points(X))
    @test maximum(abs.(r.φ .- X.φ)) < 1e-10 && abs(r.L - X.L) < 1e-12
    # a known state through another parametrization of its curve (t ↦ t + 0.3 sin t keeps 0 and 2π): the same state
    X0 = smooth_state(128); z = curve_function(X0)
    r = from_curve(t -> z(t + 0.3sin(t)), 128; tol=1e-10)
    @printf("from_curve, reparametrized: φ %.1e, L %.1e\n", maximum(abs.(r.φ .- X0.φ)), abs(r.L - X0.L))
    @test maximum(abs.(r.φ .- X0.φ)) < 1e-10 && abs(r.L - X0.L) < 1e-12 && r.x0 ≈ X0.x0
    # automatic N: the first of 64, 128, … whose φ tail is below tol/100, capped by Nmax
    a = from_curve(t -> (1 + 0.3cos(5t)) * Vec2(cos(t), sin(t)); tol=1e-6); n = length(a.φ)
    @test n > 64 && spectral_tail(a) < 1e-8 && spectral_tail(from_curve(t -> (1 + 0.3cos(5t)) * Vec2(cos(t), sin(t)), n ÷ 2)) >= 1e-8
    @test length(from_curve(γ; tol=1e-6, Nmax=64).φ) == 64
end

@testset "geometry of a state" begin
    # circle of radius R about c: A = πR², M = πR² c
    R, c = 1.3, Vec2(0.4, -0.2)
    A, M = curve_moments(ThetaLState(zeros(32), 2π * R, c + Vec2(0, -R)))
    @test A ≈ π * R^2 rtol = 1e-14
    @test M ≈ π * R^2 * c rtol = 1e-14
    circ = ThetaLState(zeros(32), 2π * R, c + Vec2(0, -R))
    @test area(circ) ≈ π * R^2 rtol = 1e-14
    @test centroid(circ) ≈ c rtol = 1e-14
    @test perimeter(circ) == 2π * R
    # a smooth curve against the shoelace formula on its trigonometric interpolant (error O(1/n²))
    X = smooth_state(64); z = curve_function(X); n = 20000
    x = [z(2π * j / n) for j in 0:n-1]; w = [x[j][1] * x[mod1(j + 1, n)][2] - x[mod1(j + 1, n)][1] * x[j][2] for j in 1:n]
    A, M = curve_moments(X)
    @test abs(A - sum(w) / 2) < 1e-7 * A
    @test norm(M - sum((x[j] + x[mod1(j + 1, n)]) * w[j] for j in 1:n) / 6) < 1e-7 * abs(A)
    # curvature: 1/R on a circle; on the ellipse 1.5 × 0.75, ab/(a² sin²t + b² cos²t)^{3/2} at the parameter t of
    # each point
    @test maximum(abs.(curvature(ThetaLState(zeros(32), 2π * R, c + Vec2(0, -R))) .- 1 / R)) < 1e-14
    e = from_curve(t -> Vec2(1.5cos(t), 0.75sin(t)), 256)
    te = [atan(x[2] / 0.75, x[1] / 1.5) for x in curve_points(e)]
    @test maximum(abs.(curvature(e) .- 1.5 * 0.75 ./ (1.5^2 .* sin.(te) .^ 2 .+ 0.75^2 .* cos.(te) .^ 2) .^ 1.5)) < 1e-9
    # arclength derivative: x_s = cos θ, and κ_s from the parametric curvature of the ellipse, κ_t/|x_t|
    @test maximum(abs.(arclength_derivative(e, first.(curve_points(e))) .- cos.(tangent_angles(e)))) < 1e-10
    κt = [-3 * 1.5 * 0.75 * (1.5^2 - 0.75^2) * sin(t) * cos(t) / (1.5^2 * sin(t)^2 + 0.75^2 * cos(t)^2)^2.5
          for t in te]
    @test maximum(abs.(arclength_derivative(e, curvature(e)) .- κt ./ hypot.(1.5 .* sin.(te), 0.75 .* cos.(te)))) < 1e-8
    # normal and tangential components: the normal itself, and the round trip
    θ = tangent_angles(X)
    @test θ ≈ 2π .* (0:63) ./ 64 .+ X.φ
    U, W = velocity_components(X, [Vec2(sin(q), -cos(q)) for q in θ])
    @test maximum(abs.(U .- 1)) < 1e-15 && maximum(abs, W) < 1e-15
    u = [Vec2(sin(3q), q^2) for q in θ]
    @test maximum(norm.(velocity_vectors(X, velocity_components(X, u)...) .- u)) < 1e-14
end

@testset "resolution: spectral tail, doubling N" begin
    X = smooth_state(64); X2 = double_grid(X)
    @test length(X2.φ) == 128 && X2.L == X.L && X2.x0 == X.x0
    @test maximum(norm.(curve_points(X2)[1:2:end] .- curve_points(X))) < 1e-13        # same curve
    g(α) = cos(3α) + 0.5cos(32α)                                                         # with the Nyquist mode
    @test maximum(abs.(double_grid(g.(2π .* (0:63) ./ 64)) .- g.(2π .* (0:127) ./ 128))) < 1e-13
    @test spectral_tail(X) < 1e-15
    @test spectral_tail(ThetaLState(X.φ .+ 1e-6 .* cos.(30 .* 2π .* (0:63) ./ 64), X.L, X.x0)) ≈ 5e-7 rtol = 1e-6
end

@testset "resolution: halving N" begin
    X = smooth_state(128); X2 = halve_grid(X)                                            # band-limited: exact
    @test length(X2.φ) == 64 && X2.L == X.L && X2.x0 == X.x0
    @test maximum(abs.(X2.φ .- X.φ[1:2:end])) < 1e-14
    @test halve_grid(double_grid(X2)).φ ≈ X2.φ atol = 1e-14
    g(α) = cos(3α) + 0.5cos(16α)                                                         # ±16 folds into the Nyquist mode of 32
    @test maximum(abs.(halve_grid(g.(2π .* (0:63) ./ 64)) .- g.(2π .* (0:31) ./ 32))) < 1e-14
end

@testset "measures: distances, integrals, the parallel curve, area rates" begin
    X1 = ThetaLState(zeros(64), 2π, Vec2(0, -1)); X2 = ThetaLState(zeros(64), π, Vec2(3, -0.5))  # unit circle; r = 1/2 at (3, 0)
    @test min_distance([X1, X2]; within=false) ≈ 1.5 rtol = 1e-12
    @test min_distance([X1]) ≈ 2sin(2π * 17 / 64 / 2) rtol = 1e-12                       # nearest points > N/4 apart
    @test curve_integral(X1, fill(1.0, 64)) ≈ 2π && abs(curve_integral(X1, first.(curve_points(X1)))) < 1e-14
    e = from_curve(t -> Vec2(1.5cos(t), 0.75sin(t)), 256); o = offset_curve(e, 0.2)
    @test o.L ≈ e.L + 2π * 0.2 rtol = 1e-12
    @test area(o) ≈ area(e) + 0.2 * e.L + π * 0.2^2 rtol = 1e-10                         # Steiner's formula
    @test o.x0 ≈ e.x0 + Vec2(0.2, 0) rtol = 1e-12                                        # n(0) = (1, 0) at (1.5, 0)
    U = [fill(2.0, 64), fill(-1.0, 64)]
    @test all(maximum(abs, u) < 1e-15 for u in remove_area_rate([X1, X2], U))
    Up = remove_area_rate([X1, X2], U; per_curve=false)                                  # Σ ∫U ds = 0 only
    @test abs(curve_integral(X1, Up[1]) + curve_integral(X2, Up[2])) < 1e-13 && Up[1][1] - Up[2][1] ≈ 3
end
