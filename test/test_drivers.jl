# The drivers (src/drivers/): `evolve` and `solve` — curve shortening of a circle (second order), the time of each
#   evaluation (U = t, exact), saving and the curves at a saved time, CSV, two curves, doubling and halving N, the
#   callback (stop, edits), a rejected step is an error; the symbol as σ(q) or σ(q, X); `dotplot` of the curves
#   (the DotPlots extension); adaptive steps (the error follows errtol, an initial layer, saved times, rejections).
#   julia -t 8 --project -e 'using Pkg; Pkg.test(test_args=["drivers"])'   (from the package root)
using ThetaL, DotPlots, LinearAlgebra, Printf, Test
using ThetaL: Vec2

circle = θ -> (cos(θ), sin(θ))

@testset "curve shortening of a circle: second order" begin
    # U = −κ with its symbol σ = −q²: R² = 1 − 2t
    e = map((20, 40)) do n
        sol = evolve(circle, (X, t) -> -curvature(X), (0, 0.45); dt=0.45 / n, symbol=q -> -q^2)
        @test sol.steps == n && sol.t == [0, 0.45]
        abs(sol.curves[end].L / 2π / sqrt(0.1) - 1)
    end
    @printf("evolve, curve shortening of a circle: R(0.45) %.2e, %.2e (ratio %.2f)\n", e..., e[1] / e[2])
    @test e[2] < 2e-3 && 3 < e[1] / e[2] < 4.5                     # 1.3e-3, ratio 3.72
end

@testset "the time of each evaluation, saving, CSV" begin
    # U = t: R = 1 + t²/2, which the midpoint rule integrates exactly if each evaluation gets its own time
    sol = evolve(circle, (X, t) -> fill(t, length(X.φ)), (0, 1); dt=0.1, saveat=0.25)
    @test sol.t == [0, 0.25, 0.5, 0.75, 1] && sol.steps == 12            # 0.1, 0.1, 0.05 to each saved time
    @test abs(sol.curves[end].L / 2π - 1.5) < 1e-14
    @test sol(0.5) === sol.curves[3] && sol(0.1 + 0.4) === sol.curves[3]
    @test_throws ArgumentError sol(0.3)
    io = IOBuffer(); save_curves(io, sol; times=[0, 1]); lines = split(String(take!(io)), "\n"; keepempty=false)
    @test lines[1] == "t,curve,x,y" && length(lines) == 1 + length(sol(0).φ) + length(sol(1).φ)
    @test parse.(Float64, split(lines[end], ",")) ≈ [1, 1, curve_points(sol(1))[end]...]
end

@testset "two curves; solve on a problem" begin
    sol = evolve([circle, θ -> (3 + 0.5cos(θ), 0.5sin(θ))],
                 (Xs, t) -> [fill(1.0, length(Xs[1].φ)), fill(-0.5, length(Xs[2].φ))], (0, 0.5); dt=0.05)
    @test abs(sol.curves[end][1].L / 2π - 1.5) < 1e-14 && abs(sol.curves[end][2].L / 2π - 0.25) < 1e-14
    @test sol.ids == fill([1, 2], 2)
    prob = CurveProblem([circle], (Xs, ev; warm) -> [-curvature(X) for X in Xs]; symbol=q -> -q^2)
    a = solve(prob, ETDRK2(); tspan=(0, 0.1), dt=0.01, saveat=0.05)
    b = evolve(circle, (X, t) -> -curvature(X), (0, 0.1); dt=0.01, saveat=0.05, symbol=q -> -q^2)
    @test a.t == b.t && a.curves[end][1].φ == b.curves[end].φ && a.curves[end][1].L == b.curves[end].L
end

@testset "doubling and halving N" begin
    # N is doubled when the curve needs it: a star given on 64 points, at rest
    star = from_curve(θ -> (1 + 0.3cos(5θ)) .* (cos(θ), sin(θ)), 64)
    sol = evolve(star, (X, t) -> zeros(length(X.φ)), (0, 0.3); dt=0.1, tol=1e-6)
    @test length(sol.curves[end].φ) > 64 && spectral_tail(sol.curves[end]) < spectral_tail(star)
    # and halved, with `coarsen`, once it has been smooth enough on half the points for 10 steps: a circle on 512
    big = ThetaLState(zeros(512), 2π, Vec2(0, -1))
    sol = evolve(big, (X, t) -> fill(0.1, length(X.φ)), (0, 0.25); dt=0.01, coarsen=true, saveat=0.05)
    @test [length(X.φ) for X in sol.curves] == [512, 512, 256, 256, 128, 128]
    @test abs(sol.curves[end].L / 2π - 1.025) < 1e-13
end

@testset "the callback: stop, edits; a rejected step is an error" begin
    stopped = evolve(circle, (X, t) -> -curvature(X), (0, 0.4); dt=0.01, saveat=0.1, symbol=q -> -q^2,
                     callback=integ -> integ.t >= 0.145 ? :stop : nothing)
    @test stopped.t[end] ≈ 0.15 && length(stopped.t) == 3
    # a callback that removes the second curve once it is small
    sol = evolve([circle, θ -> (3 + 0.3cos(θ), 0.3sin(θ))], (Xs, t) -> [-curvature(X) for X in Xs], (0, 0.1);
                 dt=0.005, saveat=0.05, symbol=q -> -q^2,
                 callback=integ -> (length(integ.states) == 2 && integ.states[2].L < 2π * 0.2 && remove_curve!(integ, 2); nothing))
    @test sol.ids[end] == [1] && length(sol.curves[end]) == 1
    @test_throws ErrorException evolve(θ -> 0.1 .* (cos(θ), sin(θ)), (X, t) -> fill(-1.0, length(X.φ)), (0, 0.3); dt=0.3)
end

@testset "the symbol: σ(q) or σ(q, X)" begin
    # the same σ given both ways, on a star (so that φ moves): the same steps
    U = (X, t) -> -curvature(X); star = θ -> (1 + 0.1cos(3θ)) .* (cos(θ), sin(θ))
    a = evolve(star, U, (0, 0.05); dt=0.01, symbol=q -> -q^2)
    b = evolve(star, U, (0, 0.05); dt=0.01, symbol=(q, X) -> -q^2 * X.L / perimeter(X))
    @test a.curves[end].φ == b.curves[end].φ && a.curves[end].L == b.curves[end].L
    @test maximum(abs, a.curves[end].φ) > 1e-3
    @test_throws ArgumentError evolve(circle, U, (0, 0.1); dt=0.02, symbol=(q, X, t) -> -q^2)
    @test_throws ArgumentError evolve(circle, U, (0, 0.1); dt=0.02, symbol=q -> -sqrt(q - 1 + 0im))
    @test_logs (:warn, r"positive rates") evolve(circle, U, (0, 0.02); dt=0.02, symbol=q -> 0.1q)
end

@testset "dotplot of curves and of an Evolution" begin
    sol = evolve(circle, (X, t) -> -curvature(X), (0, 0.1); dt=0.05, saveat=0.05, symbol=q -> -q^2)
    plot(x...; kw...) = (io = IOBuffer(); dotplot(io, x...; kw...); String(take!(io)))
    @test plot(sol) == plot([[curve_points(X)] for X in sol.curves]; t=sol.t) && !isempty(plot(sol))
    @test plot(sol; width=30) == plot([[curve_points(X)] for X in sol.curves]; t=sol.t, width=30)
    @test plot(sol(0)) == plot([[curve_points(sol(0))]])
    two = evolve([circle, θ -> (3 + 0.5cos(θ), 0.5sin(θ))], (Xs, t) -> [-curvature(X) for X in Xs], (0, 0.1);
                 dt=0.05, saveat=0.05, symbol=q -> -q^2)
    @test plot(two) == plot([curve_points.(Xs) for Xs in two.curves]; t=two.t)
    @test plot(two(0.1)) == plot([curve_points.(two(0.1))])
end

@testset "adaptive steps: the error follows errtol" begin
    # mean curvature flow of a star against a fine fixed step: the step from the estimate (err ∝ h²) makes the
    # global error of the second-order step ∝ errtol, with ~ errtol^(-1/2) steps
    star = θ -> (1 + 0.2cos(5θ)) .* (cos(θ), sin(θ)); U = (X, t) -> -curvature(X); T = 0.05
    ref = evolve(star, U, (0, T); dt=T / 4000, symbol=q -> -q^2)
    e = map((1e-4, 1e-6)) do errtol
        sol = evolve(star, U, (0, T); adaptive=true, errtol, symbol=q -> -q^2)
        @test sol.t == [0, T] && sol.rejections <= 5
        X = sol.curves[end]; Xr = ref.curves[end]
        (max(maximum(abs, X.φ .- Xr.φ), abs(X.L - Xr.L) / Xr.L, norm(X.x0 - Xr.x0) / Xr.L), sol.steps)
    end
    @printf("evolve, adaptive: error %.2e in %d steps (errtol 1e-4), %.2e in %d steps (1e-6)\n", e[1]..., e[2]...)
    @test e[1][1] < 1e-4 && e[2][1] < 1e-6                         # 1.3e-5, 1.3e-7
    @test 30 < e[1][1] / e[2][1] < 300 && 5 < e[2][2] / e[1][2] < 20 # 99; 584/60
    # no estimate, no control; no dt, no run
    @test_throws ArgumentError evolve(star, U, (0, T); adaptive=true, scheme=ExpEuler(), symbol=q -> -q^2)
    @test_throws ArgumentError evolve(star, U, (0, T); symbol=q -> -q^2)
end

@testset "adaptive steps: an initial layer, saved times, rejections" begin
    # surface diffusion of a flower (σ = −q⁴): the lobes decay within the first steps, which must be small
    flower = θ -> (1 + 0.25cos(4θ)) .* (cos(θ), sin(θ))
    κss(X) = arclength_derivative(X, arclength_derivative(X, curvature(X)))
    hs = Float64[]; ts = Float64[]; saves = [0.001, 0.0033, 0.01, 0.02]
    sol = evolve(flower, (X, t) -> κss(X), (0, 0.02); adaptive=true, errtol=1e-4, symbol=q -> -q^4, saveat=saves,
                 callback=integ -> (push!(hs, integ.stats.last_step); push!(ts, integ.t); nothing))
    @test sol.t == [0; saves]                                       # landed on exactly
    @test sol.steps == length(hs) && sol.rejections >= 1
    @test maximum(hs[1:5]) < 1e-5 && maximum(hs[ts .> 0.01]) > 100 * maximum(hs[1:5])   # ~1e-6, then ~1e-3
    @test abs(area(sol.curves[end]) - area(sol.curves[1])) < 1e-3 * area(sol.curves[1])
    # a step function bounds the step; dtmax too
    hb = Float64[]
    evolve(flower, (X, t) -> κss(X), (0, 0.005); adaptive=true, errtol=1e-4, symbol=q -> -q^4, dt=(X, t) -> 1e-4,
           dtmax=2e-4, callback=integ -> (push!(hb, integ.stats.last_step); nothing))
    @test maximum(hb) <= 1e-4 * (1 + 1e-12) && maximum(hb) > 5e-5
    # a velocity that throws StepFailure for steps over 0.01 (its own reason): the step is shrunk and retried
    reasons = Any[]
    prob = CurveProblem([circle], (Xs, ev; warm) -> ev.h > 0.01 ? throw(StepFailure(:too_long)) :
                                                  [fill(0.1, length(X.φ)) for X in Xs])
    hf = Float64[]
    sol = solve(prob; tspan=(0, 0.5), adaptive=true, dt=0.1, saveat=0.25,
                callback=integ -> (push!(hf, integ.stats.last_step); push!(reasons, integ.stats.last_rejection); nothing))
    @test sol.t == [0, 0.25, 0.5] && sol.rejections > 0 && maximum(hf) <= 0.01 && :too_long in reasons
    @test abs(sol.curves[end][1].L / 2π - 1.05) < 1e-13
    # always rejected: an error, not a loop
    never = CurveProblem([circle], (Xs, ev; warm) -> ev.kind === :init ? [zeros(length(X.φ)) for X in Xs] :
                                                     throw(StepFailure(:never)))
    @test_throws ErrorException solve(never; tspan=(0, 1), adaptive=true)
    @test_throws ErrorException solve(never; tspan=(0, 1), adaptive=true, dtmin=1e-3)
end
