# The integrator (src/integrator/, src/problem/): the velocity protocol (the sequence of evaluations, their kinds,
#   stages, times and contexts; results as a vector or a named tuple; wrong sizes; a plain tuple (an error)), exact
#   landing, curve edits under a coupled velocity (U invalidated, ids stable; a removal equals a fresh start),
#   rejection (StepFailure, a non-positive length, the error estimate before the end evaluation, `validate`) and
#   the retry's context, a scheme change, the automatic filter cutoff, checkpoints (round trip), aliasing, and a
#   coupled deforming pair against an accurate ODE.
#   julia -t 8 --project -e 'using Pkg; Pkg.test(test_args=["integrator"])'   (from the package root)
using ThetaL, LinearAlgebra, Printf, Test
using ThetaL: Vec2
include(joinpath(@__DIR__, "common", "curves.jl"))

circle(c, r) = θ -> (c[1] + r * cos(θ), c[2] + r * sin(θ))
radius(X) = X.L / 2π
# particles coupled through a mean field: V_i = κ̄ − κ_i, κ̄ = 2π·(number of curves)/ΣL (total area kept)
mean_field(states) = (κbar = 2π * length(states) / sum(X -> X.L, states); [κbar .- curvature(X) for X in states])

@testset "the evaluation sequence: kinds, stages, times, contexts" begin
    log = []
    function velocity(states, ev; warm)
        push!(log, (ev.kind, ev.stage, ev.t, warm, ev.attempt, copy(ev.ids)))
        (; U=[-curvature(X) for X in states], ctx=length(log))
    end
    integ = init(CurveProblem([circle((0, 0), 1)], velocity; symbol=q -> -q^2), ETDRK2(); t0=0.5)
    @test advance!(integ, 0.6) === :success && advance!(integ, 0.7) === :success
    @test [(l[1], l[2], l[4]) for l in log] == [(:init, 0, nothing), (:stage, 2, 1), (:end, 3, 2), (:stage, 2, 3), (:end, 3, 4)]
    @test log[1][3] == 0.5 && log[2][3] ≈ 0.55 && log[3][3] === 0.6 && log[5][3] === 0.7     # exact landing
    @test integ.t === 0.7 && integ.ctx == 5 && integ.stats.evaluations == 5 && integ.stats.steps == 2
    # ExpEuler: one evaluation per step
    integ = init(CurveProblem([circle((0, 0), 1)], velocity; symbol=q -> -q^2), ExpEuler())
    n = length(log); step!(integ, 0.1); step!(integ, 0.1)
    @test [l[1] for l in log[n+1:end]] == [:end, :end] && integ.t ≈ 0.2
    @test_throws ArgumentError advance!(integ, integ.t)                      # not after the current time
end

@testset "velocity results: vector, named tuple, wrong sizes, a plain tuple (an error)" begin
    X = from_curve(circle((0, 0), 1), 64)
    ok(ret) = init(CurveProblem([X], (Xs, ev; warm) -> ret(Xs)), ETDRK2())
    @test ok(Xs -> [zeros(64)]).ctx === nothing
    @test ok(Xs -> (; U=[zeros(64)], ctx=:mine, ydot=1.0)).ctx === :mine       # unknown fields are accepted
    @test_throws DimensionMismatch ok(Xs -> [zeros(63)])
    @test_throws DimensionMismatch ok(Xs -> [zeros(64), zeros(64)])
    @test_throws ArgumentError ok(Xs -> ([zeros(64)], [zeros(64)]))
    @test_throws ArgumentError evolve(X, (X, t) -> (zeros(64), zeros(64)), (0, 0.1); dt=0.1)
end

@testset "curve edits: velocities invalidated, ids stable, a removal equals a fresh start" begin
    curves = [circle((0, 0), 1.0), circle((4, 0), 0.6), circle((0, 4), 0.8)]
    prob = CurveProblem(curves, (Xs, ev; warm) -> mean_field(Xs); symbol=q -> -q^2, tags=[:a, :b, :c])
    integ = init(prob, ETDRK2()); advance!(integ, 0.01)
    after_one = copy(integ.states)
    remove_curve!(integ, 2)
    @test integ.U === nothing && integ.ids == [1, 3] && integ.tags == [:a, :c] && integ.topology == 1
    advance!(integ, 0.02)
    fresh = init(CurveProblem(after_one[[1, 3]], (Xs, ev; warm) -> mean_field(Xs); symbol=q -> -q^2), ETDRK2();
                 t0=0.01)
    advance!(fresh, 0.02)
    @test all(a.φ == b.φ && a.L == b.L && a.x0 == b.x0 for (a, b) in zip(integ.states, fresh.states))
    # add, set, offset: each invalidates; new ids continue; double and halve keep U
    id = add_curve!(integ, circle((4, 4), 0.5); tag=:d)
    @test id == 4 && integ.ids == [1, 3, 4] && integ.U === nothing && curve_index(integ, 4) == 3
    advance!(integ, 0.03); n = length(integ.states[1].φ)
    double_grid!(integ, 1)
    @test integ.U !== nothing && length(integ.U[1]) == 2n && length(integ.states[1].φ) == 2n
    halve_grid!(integ, 1)
    @test integ.U !== nothing && length(integ.U[1]) == n
    offset_curve!(integ, 2, 0.1)
    @test integ.U === nothing && integ.topology == 5
    set_curve!(integ, 3, circle((4, 4), 0.4))
    @test integ.U === nothing && advance!(integ, 0.04) === :success
end

@testset "rejection: StepFailure, a non-positive length, the error estimate, validate; the retry's context" begin
    # a solver that fails once at the midpoint of its second step
    log = []
    function velocity(states, ev; warm)
        push!(log, (ev.kind, ev.attempt, warm))
        ev.kind === :stage && ev.t > 0.1 && ev.attempt == 1 && throw(StepFailure(:not_converged))
        (; U=[-curvature(X) for X in states], ctx=length(log))
    end
    integ = init(CurveProblem([circle((0, 0), 1)], velocity; symbol=q -> -q^2), ETDRK2())
    advance!(integ, 0.1)
    before = (integ.t, integ.states, integ.U, integ.ctx)
    @test advance!(integ, 0.2) === :rejected
    @test (integ.t, integ.states, integ.U, integ.ctx) === before && integ.stats.last_rejection === :not_converged
    @test integ.stats.rejections == 1 && integ.stats.consecutive_rejections == 1
    @test advance!(integ, 0.2) === :success
    @test log[end-1] == (:stage, 2, before[4])                          # the retry starts from the committed context
    @test integ.stats.consecutive_rejections == 0 && integ.t === 0.2
    # a shrinking circle whose midpoint stage would have no length left
    shrink = init(CurveProblem([circle((0, 0), 0.1)], (Xs, ev; warm) -> [fill(-1.0, length(X.φ)) for X in Xs]), ETDRK2())
    @test advance!(shrink, 0.3) === :rejected && shrink.stats.last_rejection === :nonpositive_length && shrink.t == 0
    # the error estimate, checked before the end evaluation: one evaluation spent. (Surface diffusion of a flower:
    # in its initial layer the estimate is first order, 2e-2 at Δt = 1e-4, 1e-7 at 1e-8.)
    κss(X) = arclength_derivative(X, arclength_derivative(X, curvature(X)))
    flower = θ -> (1 + 0.25cos(4θ)) .* (cos(θ), sin(θ))
    sd = init(CurveProblem([flower], (Xs, ev; warm) -> κss.(Xs); symbol=q -> -q^4), ETDRK2())
    n = sd.stats.evaluations
    @test advance!(sd, 1e-4; errtol=1e-3) === :rejected
    @test sd.stats.last_rejection === :error_estimate && sd.stats.evaluations == n + 1 && sd.stats.error_estimate > 1e-3
    @test advance!(sd, 1e-8; errtol=1e-6) === :success
    # validate
    prob = CurveProblem([circle((0, 0), 1)], (Xs, ev; warm) -> [fill(-1.0, length(X.φ)) for X in Xs];
                        validate=Xs -> all(X -> X.L > 2π * 0.95, Xs))
    v = init(prob, ETDRK2())
    @test advance!(v, 0.01) === :success && advance!(v, 0.1) === :rejected && v.stats.last_rejection === :invalid_state
end

@testset "scheme changes, the automatic filter cutoff" begin
    prob = CurveProblem([circle((0, 0), 1)], (Xs, ev; warm) -> [-curvature(X) for X in Xs]; symbol=q -> -q^2)
    integ = init(prob, ETDRK2()); advance!(integ, 0.01); n = integ.stats.evaluations
    integ.scheme = ExpEuler(); advance!(integ, 0.02)
    @test integ.stats.evaluations == n + 1
    auto = init(prob, ETDRK2(filter_cutoff=:auto))
    @test auto.scheme.filter.cutoff isa Int && auto.scheme.filter.cutoff >= 32
end

@testset "checkpoints: plain data, a round trip that continues identically" begin
    curves = [circle((0, 0), 1.0), θ -> (3 + 0.5cos(θ) + 0.1cos(3θ), 0.5sin(θ))]
    prob = CurveProblem(curves, (Xs, ev; warm) -> mean_field(Xs); symbol=q -> -q^2, tags=[:outer, "second"])
    a = init(prob, ETDRK2(filter_cutoff=:auto))
    for t in (0.01, 0.02, 0.03); advance!(a, t); end
    remove_curve!(a, 1); add_curve!(a, circle((0, 0), 0.9); tag=:new); advance!(a, 0.04)
    cp = checkpoint(a)
    @test cp.format == 1 && cp.phi isa Vector{Vector{Float64}} && cp.x0[1] isa Vector{Float64}
    @test cp.ids == [2, 3] && cp.tags == ["second", :new] && cp.scheme.filter_cutoff isa Int
    b = init(prob, ETDRK2(filter_cutoff=:auto); from=cp)
    @test b.t == a.t && b.ids == a.ids && b.topology == a.topology && b.scheme.filter == a.scheme.filter
    for t in (0.05, 0.06); advance!(a, t); advance!(b, t); end
    @test all(x.φ == y.φ && x.L == y.L && x.x0 == y.x0 for (x, y) in zip(a.states, b.states))
    @test add_curve!(a, circle((5, 5), 0.3)) == add_curve!(b, circle((5, 5), 0.3)) == 4
end

@testset "aliasing: what a caller keeps stays valid" begin
    integ = init(CurveProblem([circle((0, 0), 1)], (Xs, ev; warm) -> [-curvature(X) for X in Xs]; symbol=q -> -q^2))
    states, U = integ.states, integ.U; φ = copy(states[1].φ); u = copy(U[1])
    advance!(integ, 0.01); double_grid!(integ, 1); advance!(integ, 0.02)
    @test states[1].φ == φ && U[1] == u && length(states[1].φ) == length(φ)
    # a velocity that reuses its buffer cannot change the integrator's data
    buf = [zeros(64)]
    reuse = init(CurveProblem([from_curve(circle((0, 0), 1), 64)], (Xs, ev; warm) -> (buf[1] .= -1.0; buf)))
    kept = reuse.U[1]; buf[1] .= 5.0
    @test all(kept .== -1.0)
end

@testset "a coupled deforming pair against an accurate ODE" begin
    # two circles under V_i = κ̄ − κ_i stay circles; their radii follow R_i′ = 2/(R₁ + R₂) − 1/R_i (area kept)
    f(R) = (m = 2 / (R[1] + R[2]); [m - 1 / R[1], m - 1 / R[2]])
    function ode(R, T, n)
        h = T / n
        for _ in 1:n
            k1 = f(R); k2 = f(R .+ h / 2 .* k1); k3 = f(R .+ h / 2 .* k2); k4 = f(R .+ h .* k3)
            R = R .+ h / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        end
        R
    end
    T = 0.1; Rex = ode([1.0, 0.6], T, 10_000)
    prob = CurveProblem([circle((0, 0), 1.0), circle((4, 0), 0.6)], (Xs, ev; warm) -> mean_field(Xs);
                        symbol=q -> -q^2)
    e = map((0.01, 0.005)) do dt
        sol = solve(prob, ETDRK2(); tspan=(0, T), dt)
        maximum(abs.(radius.(sol.curves[end]) .- Rex))
    end
    @printf("coupled pair: radius errors %.2e %.2e, ratio %.2f\n", e..., e[1] / e[2])
    @test e[2] < 1e-5 && 3.5 < e[1] / e[2] < 4.5
end
