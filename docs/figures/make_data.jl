# The data behind the figures of ThetaL's documentation, as CSV files in figures/data/ (plot_figures.py reads
# them). Every run is small (a few minutes in all).
#   julia -t 8 --project=docs docs/figures/make_data.jl          (from the package root)
using ThetaL, LinearAlgebra, Printf
using ThetaL: Vec2, Lawson2

const DATA = joinpath(@__DIR__, "data")
mkpath(DATA)
csv(name, header, rows) = open(joinpath(DATA, name), "w") do io
    println(io, join(header, ","))
    for r in rows
        println(io, join(r, ","))
    end
    println("wrote data/", name, " (", length(rows), " rows)")
end
# the points of a curve for a figure: at most 512 (every k-th grid point), the first repeated to close it
points_rows(t, id, X) = (x = curve_points(X); k = max(1, length(x) ÷ 512); [(t, id, p[1], p[2]) for p in [x[1:k:end]; [x[1]]]])
field(u) = (X, t) -> velocity_components(X, u.(curve_points(X)))[1]
κss(X) = arclength_derivative(X, arclength_derivative(X, curvature(X)))

# ---- gallery: curves at saved times -----------------------------------------------------------------------------
star = θ -> (1 + 0.3cos(5θ)) .* (cos(θ), sin(θ))
T = area(from_curve(star)) / 2π
sol = evolve(star, (X, t) -> -curvature(X), (0, 0.95T); dt=(X, t) -> 0.0125 / maximum(abs2, curvature(X)),
             symbol=q -> -q^2, saveat=0.95T / 8)
csv("gallery_mcf.csv", ("t", "curve", "x", "y"), reduce(vcat, [points_rows(t, 1, X) for (t, X) in zip(sol.t, sol.curves)]))
csv("mcf_area.csv", ("t", "area", "exact"), [(t, area(X), area(sol.curves[1]) - 2π * t) for (t, X) in zip(sol.t, sol.curves)])

flower = θ -> (1 + 0.25cos(4θ)) .* (cos(θ), sin(θ))
sol = evolve(flower, (X, t) -> κss(X), (0, 0.05); dt=2.5e-5, symbol=q -> -q^4, saveat=[0, 0.001, 0.003, 0.01, 0.02, 0.05])
csv("gallery_sd.csv", ("t", "curve", "x", "y"), reduce(vcat, [points_rows(t, 1, X) for (t, X) in zip(sol.t, sol.curves)]))

ω(r) = r < 1e-8 ? 1 - r^2 / 2 : -expm1(-r^2) / r^2
sol = evolve(θ -> (1.0 + 0.5cos(θ), 0.5sin(θ)), field(x -> ω(hypot(x...)) .* (-x[2], x[1])), (0, 10.0);
             dt=(X, t) -> 0.5perimeter(X) / length(X.φ), saveat=[0, 2.5, 5, 10])
csv("gallery_vortex.csv", ("t", "curve", "x", "y"), reduce(vcat, [points_rows(t, 1, X) for (t, X) in zip(sol.t, sol.curves)]))
csv("vortex_resolution.csv", ("t", "N", "length"), [(t, length(X.φ), X.L) for (t, X) in zip(sol.t, sol.curves)])

stiffness(ϑ) = 1 - 15 * 0.05 * cos(4ϑ)
nangles(X) = tangent_angles(X) .- π / 2
prob = CurveProblem([θ -> (cos(θ), sin(θ))], (Xs, ev; warm) -> [-stiffness.(nangles(X)) .* curvature(X) for X in Xs];
                    symbol=(X, k, info) -> maximum(stiffness, nangles(X)) .* (2π / perimeter(X) .* abs.(k)) .^ 2)
sol = solve(prob, ETDRK2(); tspan=(0, 0.45), dt=(Xs, t) -> 0.01 / maximum(abs2, curvature(Xs[1])), saveat=0.05)
csv("gallery_aniso.csv", ("t", "curve", "x", "y"), reduce(vcat, [points_rows(t, 1, X[1]) for (t, X) in zip(sol.t, sol.curves)]))
csv("aniso_area.csv", ("t", "area", "exact"), [(t, area(X[1]), π - 2π * t) for (t, X) in zip(sol.t, sol.curves)])

# coarsening: twelve particles under a mean field, removed as they vanish
particles = [((3.0i, 3.0j), 0.3 + 0.7 * mod(0.618i + 0.382j * i + 0.1j, 1.0)) for i in 1:4 for j in 1:3]
mf(Xs, ev; warm) = (κbar = 2π * length(Xs) / sum(perimeter, Xs);
                    remove_area_rate(Xs, [κbar .- curvature(X) for X in Xs]; per_curve=false))
integ = init(CurveProblem([θ -> (c[1] + r * cos(θ), c[2] + r * sin(θ)) for (c, r) in particles], mf; symbol=q -> -q^2))
shapes = []; radii = []; tnext = 0.0
function record!(integ)
    for (id, X) in zip(integ.ids, integ.states)
        push!(radii, (integ.t, id, perimeter(X) / 2π))
    end
    if integ.t >= tnext
        for (id, X) in zip(integ.ids, integ.states); append!(shapes, points_rows(integ.t, id, X)); end
        global tnext += 1.0
    end
end
record!(integ)
while integ.t < 3 && length(integ.states) > 1
    h = 0.05 * minimum(X -> (perimeter(X) / 2π)^2, integ.states)
    advance!(integ, integ.t + h) === :success || (remove_curve!(integ, argmin(perimeter.(integ.states))); continue)
    for i in reverse(eachindex(integ.states))
        perimeter(integ.states[i]) / 2π < 0.05 && remove_curve!(integ, i)
    end
    record!(integ)
end
csv("gallery_coarsening.csv", ("t", "curve", "x", "y"), shapes)
csv("coarsening_radii.csv", ("t", "curve", "radius"), radii)

# ---- convergence ------------------------------------------------------------------------------------------------
# rigid rotation of a smooth curve about an outside point, one turn: the rigid part taken out (any scheme) or not
X0 = ThetaLState([0.3sin(2α) + 0.1cos(4α) for α in 2π .* (0:127) ./ 128], 2π, Vec2(0.2, -1.0)); z0 = curve_function(X0)
function curve_distance(z, xs)
    αs = range(0, 2π, length=4097)[1:end-1]; zs = z.(αs); h = 1e-4
    maximum(xs) do x
        a = αs[argmin([norm(y - x) for y in zs])]
        for _ in 1:8
            d1 = (z(a + h) - z(a - h)) / 2h; d2 = (z(a + h) - 2z(a) + z(a - h)) / h^2
            a -= dot(z(a) - x, d1) / (dot(d1, d1) + dot(z(a) - x, d2))
        end
        norm(z(a) - x)
    end
end
rows = []
for dt in (0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 5e-4, 2.5e-4)
    for (name, scheme, rigid) in (("ETDRK2", ETDRK2(), :full), ("ExpEuler", ExpEuler(), :full), ("ETDRK2 without the rigid split", ETDRK2(), :none))
        e = try
            X = evolve(X0, field(x -> 2π * Vec2(-x[2] + 0.5, x[1] - 2)), (0, 1); dt, scheme, rigid, Nmax=128).curves[end]
            all(isfinite, X.φ) ? curve_distance(z0, curve_points(X)) : NaN
        catch
            NaN
        end
        push!(rows, (dt, name, isfinite(e) && e < 10 ? e : NaN))
    end
end
csv("conv_rigid.csv", ("dt", "scheme", "error"), rows)

# a linear flow u = A x: the exact curve is exp(tA) applied to the first
A = [0.3 0.5; -0.4 -0.1]; Tl = 0.4; back = x -> Vec2(exp(-Tl * A) * x)
rows = []
for dt in (0.1, 0.05, 0.025, 0.0125), (name, scheme) in (("ETDRK2", ETDRK2()), ("ExpEuler", ExpEuler()))
    X = evolve(X0, field(x -> Vec2(A * x)), (0, Tl); dt, scheme, symbol=q -> -0.5q).curves[end]
    push!(rows, (dt, name, curve_distance(z0, back.(curve_points(X)))))
end
csv("conv_linear.csv", ("dt", "scheme", "error"), rows)

# surface diffusion: from the shape itself (initial layer) and from a relaxed state, against a fine reference
sdV(X, t) = κss(X); sd = q -> -q^4
raw = from_curve(θ -> (1 + 0.1cos(3θ)) .* (cos(θ), sin(θ)), 64)
relaxed = evolve(raw, sdV, (0, 0.005); dt=1e-5, symbol=sd, Nmax=64).curves[end]
rows = []
for (start, X1) in (("from the shape", raw), ("relaxed", relaxed)), (name, scheme) in (("ETDRK2", ETDRK2()), ("ExpEuler", ExpEuler()))
    ref = evolve(X1, sdV, (0, 0.05); dt=2.5e-5, symbol=sd, Nmax=64, scheme).curves[end]
    for dt in (4e-3, 2e-3, 1e-3, 5e-4, 2.5e-4)
        X = evolve(X1, sdV, (0, 0.05); dt, symbol=sd, Nmax=64, scheme).curves[end]
        push!(rows, (dt, name, start, abs(X.L - ref.L)))
    end
end
csv("conv_sd.csv", ("dt", "scheme", "start", "error"), rows)

# a manufactured deforming ellipse under wrong symbols
ωe = 1.0; a(t) = 1.5 + 0.3sin(t); b(t) = 0.8 - 0.2sin(2t); ȧ(t) = 0.3cos(t); ḃ(t) = -0.4cos(2t)
rot(x, t) = (cs = cos(ωe * t); sn = sin(ωe * t); Vec2(cs * x[1] + sn * x[2], -sn * x[1] + cs * x[2]))
function level(x, t)
    p = rot(x, t); F = (p[1] / a(t))^2 + (p[2] / b(t))^2 - 1
    g = Vec2(2p[1] / a(t)^2, 2p[2] / b(t)^2)
    F, -2p[1]^2 * ȧ(t) / a(t)^3 - 2p[2]^2 * ḃ(t) / b(t)^3 + dot(g, Vec2(ωe * p[2], -ωe * p[1])), norm(g)
end
Vm(X, t) = [(r = level(x, t); -r[2] / r[3]) for x in curve_points(X)]
e0 = from_curve(t -> Vec2(a(0) * cos(t), b(0) * sin(t)), 256)
rows = []
for scale in (0.0, 0.5, 2.0), dt in (0.01, 0.005, 0.0025, 0.00125)
    X = evolve(e0, Vm, (0, 1.0); dt, symbol=q -> -scale * q, Nmax=256).curves[end]
    push!(rows, (dt, scale, maximum(x -> (r = level(x, 1.0); abs(r[1]) / r[3]), curve_points(X))))
end
csv("conv_manufactured.csv", ("dt", "symbol_scale", "error"), rows)

# step-size control on surface diffusion of a flower: the accepted steps and their estimates
fl = init(CurveProblem([flower], (Xs, ev; warm) -> κss.(Xs); symbol=q -> -q^4), ETDRK2())
rows = []; h = 1e-8; tol = 1e-4
while fl.t < 0.2
    t0 = fl.t
    st = advance!(fl, min(fl.t + h, 0.2); errtol=tol)
    st === :success && push!(rows, (fl.t, fl.t - t0, fl.stats.error_estimate))
    global h *= clamp(0.9 * sqrt(tol / fl.stats.error_estimate), 0.2, 5.0)
end
csv("step_control.csv", ("t", "dt", "error_estimate"), rows)
@printf("step control: %d steps, %d rejected\n", fl.stats.steps, fl.stats.rejections)
