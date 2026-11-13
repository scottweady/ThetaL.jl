# State carried with the curves (src/problem/carried.jl, src/stepping/carried.jl): a bulk state driving a circle
#   (exact with its symbol, second order without), fields on the curves under rigid motion (unchanged, every
#   transport), a conserved density under a deforming flow (∫c ds kept), a material field against its exact image
#   (second order), surface diffusion of a field treated exactly by its symbol, and the bookkeeping: rejection
#   rolls the carried state back, the curve edits carry the fields, checkpoints continue identically, `solve` keeps
#   fields and bulk, missing W or ydot is an error.
#   julia -t 8 --project -e 'using Pkg; Pkg.test(test_args=["carried"])'   (from the package root)
using ThetaL, LinearAlgebra, Printf, Test
using ThetaL: Vec2
include(joinpath(@__DIR__, "common", "curves.jl"))

grid(N) = 2π .* (0:N-1) ./ N
unit_circle(N) = ThetaLState(zeros(N), 2π, Vec2(0, -1))
# the normal and tangential velocities of a velocity field u(x) on every curve, as a velocity's named tuple
flow(u) = Xs -> (UW = [velocity_components(X, u.(curve_points(X))) for X in Xs]; (U=first.(UW), W=last.(UW)))

@testset "a bulk state driving a circle" begin
    # U = −y on the unit circle with y' = −y: y = y₀e^{−t}, R = 1 − y₀(1 − e^{−t}). With the symbol σ = −1 the bulk
    # is exact at any step; the radius (L by the step's quadrature of U) is second order
    y0, T = 0.5, 1.0
    velocity(Xs, ev; warm) = (; U=[fill(-ev_y(ev), length(X.φ)) for X in Xs], ydot=[-ev_y(ev)])
    ev_y(ev) = ev.bulk[1]
    for (σ, exact_bulk) in ((nothing, false), ([-1.0], true))
        e = map((0.1, 0.05)) do dt
            prob = CurveProblem([unit_circle(64)], velocity; bulk=BulkState([y0]; symbol=σ))
            sol = solve(prob; tspan=(0, T), dt=dt)
            ey = abs(sol.bulk[end][1] - y0 * exp(-T))
            exact_bulk && @test ey < 1e-15
            (ey, abs(sol.curves[end][1].L / 2π - (1 - y0 * (1 - exp(-T)))))
        end
        @printf("bulk-driven circle, bulk symbol %s: y errors %.1e %.1e, R errors %.2e %.2e (ratio %.2f)\n",
                σ === nothing ? "none" : "-1", e[1][1], e[2][1], e[1][2], e[2][2], e[1][2] / e[2][2])
        @test 3.5 < e[1][2] / e[2][2] < 4.5
        exact_bulk || @test 3.5 < e[1][1] / e[2][1] < 4.5
    end
end

@testset "fields under rigid motion are unchanged" begin
    # a smooth curve turning about a point outside it and drifting; the material moves with it, so W − T = 0
    u(x) = 2π * Vec2(-(x[2] + 0.5), x[1] - 2) + Vec2(0.3, -0.2)
    X0 = smooth_state(128); c0 = 1 .+ 0.3 .* cos.(3 .* grid(128))
    for transport in (:density, :material, :none)
        velocity(Xs, ev; warm) = merge(flow(u)(Xs), (; rates=(; c=[zeros(length(X.φ)) for X in Xs])))
        prob = CurveProblem([X0], velocity; fields=(; c=CurveField([c0]; transport)))
        sol = solve(prob; tspan=(0, 1), dt=0.1, tol=1.0)
        @test maximum(abs.(sol.fields[end].c[1] .- c0)) < 1e-12
    end
end

@testset "a conserved density under a deforming flow" begin
    # u = Ax deforms the curve; ∫c ds = (L/N) Σc is kept to roundoff by the conservation form
    A = [0.3 0.5; -0.4 -0.1]; u(x) = Vec2(A * x)
    X0 = smooth_state(128); c0 = 1 .+ 0.3 .* cos.(3 .* grid(128)) .+ 0.1 .* sin.(5 .* grid(128))
    prob = CurveProblem([X0], (Xs, ev; warm) -> flow(u)(Xs); fields=(; c=CurveField([c0]; transport=:density)))
    sol = solve(prob; tspan=(0, 1), dt=0.02, tol=1.0)
    mass(X, c) = X.L / length(c) * sum(c)
    drift = abs(mass(sol.curves[end][1], sol.fields[end].c[1]) / mass(X0, c0) - 1)
    @printf("density under a linear flow: ∫c ds drift %.1e, L %.3f → %.3f\n", drift, X0.L, sol.curves[end][1].L)
    @test drift < 1e-12 && abs(sol.curves[end][1].L - X0.L) > 0.1
end

@testset "a material field against its exact image" begin
    # u = Ax moves material points as x(t) = e^{At}x₀, so c = g(x₀) is g(e^{−At}x) at the points x of the curve
    A = [0.3 0.5; -0.4 -0.1]; u(x) = Vec2(A * x); g(x) = x[1] + 0.5x[2]^2; T = 1.0
    X0 = smooth_state(128)
    prob = CurveProblem([X0], (Xs, ev; warm) -> flow(u)(Xs);
                        fields=(; c=CurveField(X -> g.(curve_points(X)); transport=:material)))
    e = map((0.04, 0.02, 0.01)) do dt
        sol = solve(prob; tspan=(0, T), dt=dt, tol=1.0); B = exp(-A * T)
        maximum(abs(c - g(B * x)) for (c, x) in zip(sol.fields[end].c[1], curve_points(sol.curves[end][1])))
    end
    @printf("material field under a linear flow: errors %.2e %.2e %.2e, ratios %.2f %.2f\n", e..., e[1] / e[2],
            e[2] / e[3])
    @test 3.5 < e[1] / e[2] < 4.5 && 3.5 < e[2] / e[3] < 4.5
end

@testset "surface diffusion of a field, exact with its symbol" begin
    # c_t = D c_ss on the unit circle at rest: c = cos 3α e^{−9Dt}; with σ = −Dq² the step is exact at D·q²·Δt = 4.5
    D, N = 1.0, 64
    function velocity(Xs, ev; warm)
        S = [D .* arclength_derivative(X, arclength_derivative(X, c)) for (X, c) in zip(Xs, ev.fields.c)]
        (; U=[zeros(length(X.φ)) for X in Xs], rates=(; c=S))
    end
    prob = CurveProblem([unit_circle(N)], velocity;
                        fields=(; c=CurveField([cos.(3 .* grid(N))]; transport=:none, symbol=q -> -D * q^2)))
    sol = solve(prob; tspan=(0, 1), dt=0.5)
    err = maximum(abs.(sol.fields[end].c[1] .- exp(-9D) .* cos.(3 .* grid(N)))) / exp(-9D)
    @printf("surface diffusion of a field, Dq²Δt = 4.5: relative error %.1e\n", err)
    @test err < 1e-10
end

@testset "bookkeeping: rejection, edits, checkpoints, solve, errors" begin
    A = [0.3 0.5; -0.4 -0.1]; u(x) = Vec2(A * x); fail = Ref(false)
    function velocity(Xs, ev; warm)
        fail[] && ev.kind === :stage && ev.attempt == 1 && throw(StepFailure(:test))
        merge(flow(u)(Xs), (; rates=(; c=[0.1 .* c for c in ev.fields.c]), ydot=-ev.bulk))
    end
    prob = CurveProblem([smooth_state(64), ξ -> (3 + 0.5cos(ξ), 0.5sin(ξ))], velocity;
                        fields=(; c=CurveField(X -> 1 .+ 0.2 .* cos.(2 .* grid(length(X.φ))))),
                        bulk=BulkState([1.0, 2.0]; symbol=[-1.0, -1.0]))
    integ = init(prob)
    @test integ.fields.c[1] ≈ 1 .+ 0.2 .* cos.(2 .* grid(64)) && integ.bulk == [1.0, 2.0]
    # a rejected step leaves the carried state as it was
    f0, y0 = integ.fields, integ.bulk
    fail[] = true
    @test advance!(integ, 0.01) === :rejected && integ.fields === f0 && integ.bulk === y0
    @test advance!(integ, 0.01) === :success && integ.fields !== f0 && integ.bulk ≈ y0 .* exp(-0.01)
    # edits: doubling interpolates, adding uses the field's init, removing drops, replacing resamples
    c1 = integ.fields.c[1]; double_grid!(integ, 1)
    @test length(integ.fields.c[1]) == 128 && maximum(abs.(integ.fields.c[1][1:2:end] .- c1)) < 1e-12
    fail[] = false
    id = add_curve!(integ, ξ -> (-3 + 0.4cos(ξ), 0.4sin(ξ)))
    @test length(integ.fields.c) == 3 && integ.fields.c[3] ≈ 1 .+ 0.2 .* cos.(2 .* grid(length(integ.states[3].φ)))
    remove_curve!(integ, 2)
    @test length(integ.fields.c) == 2 && integ.ids == [1, id]
    set_curve!(integ, 1, halve_grid(integ.states[1]))
    @test length(integ.fields.c[1]) == 64 && maximum(abs.(integ.fields.c[1] .- c1)) < 1e-12
    @test advance!(integ, 0.02) === :success
    # a checkpoint continues identically, carried state included
    cp = checkpoint(integ); b = init(prob; from=cp)
    for t in (0.03, 0.04); advance!(integ, t); advance!(b, t); end
    @test all(x.φ == y.φ for (x, y) in zip(integ.states, b.states))
    @test integ.fields.c == b.fields.c && integ.bulk == b.bulk
    # solve keeps the fields and the bulk at the saved times
    sol = solve(prob; tspan=(0, 0.04), dt=0.01, saveat=0.02)
    @test length(sol.fields) == length(sol.bulk) == 3 && sol.bulk[end] ≈ [1.0, 2.0] .* exp(-0.04)
    @test length(sol.fields[end].c) == 2
    # missing W (a transported field) or ydot (a bulk state) is an error
    @test_throws ArgumentError init(CurveProblem([smooth_state(64)], (Xs, ev; warm) -> [zeros(64)];
                                                 fields=(; c=CurveField(1.0))))
    @test_throws ArgumentError init(CurveProblem([smooth_state(64)], (Xs, ev; warm) -> [zeros(64)];
                                                 bulk=BulkState([1.0])))
end
