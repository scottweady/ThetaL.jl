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

# ---- the examples: frames for the movies, and their checks -----------------------------------------------------------------------------
star = θ -> (1 + 0.3cos(5θ)) .* (cos(θ), sin(θ))
T = area(from_curve(star)) / 2π
sol = evolve(star, (X, t) -> -curvature(X), (0, 0.95T); dt=(X, t) -> 0.0125 / maximum(abs2, curvature(X)),
             symbol=q -> -q^2, saveat=0.95T / 64)
csv("mcf_movie.csv", ("t", "curve", "x", "y"), reduce(vcat, [points_rows(t, 1, X) for (t, X) in zip(sol.t, sol.curves)]))
csv("mcf_area.csv", ("t", "area", "exact"), [(t, area(X), area(sol.curves[1]) - 2π * t) for (t, X) in zip(sol.t[1:8:end], sol.curves[1:8:end])])

flower = θ -> (1 + 0.25cos(4θ)) .* (cos(θ), sin(θ))        # for the step control below

ω(r) = r < 1e-8 ? 1 - r^2 / 2 : -expm1(-r^2) / r^2
sol = evolve(θ -> (1.0 + 0.5cos(θ), 0.5sin(θ)), field(x -> ω(hypot(x...)) .* (-x[2], x[1])), (0, 10.0);
             dt=(X, t) -> 0.5perimeter(X) / length(X.φ), saveat=0:0.125:10)
csv("vortex_movie.csv", ("t", "curve", "x", "y"), reduce(vcat, [points_rows(t, 1, X) for (t, X) in zip(sol.t, sol.curves)]))
csv("vortex_resolution.csv", ("t", "N", "length"), [(t, length(X.φ), X.L) for (t, X) in zip(sol.t, sol.curves)])

# Hele-Shaw: an expanding bubble, the exterior Dirichlet problem by a Cauchy integral (examples/hele_shaw.jl)
function hele_shaw(X, t; S=0.005, J=2π)
    N, ds = length(X), perimeter(X) / length(X)
    z = [complex(x...) for x in curve_points(X)]; τ = cis.(tangent_angles(X)); zc = complex(centroid(X)...)
    C = [i == j ? 0.0im : τ[j] * ds / (2π * im * (z[j] - z[i])) for i in 1:N, j in 1:N]
    c = vec(sum(C; dims=2))
    A = real.(C); A[diagind(A)] .-= real.(c); A .+= ds
    μ = A \ (-S .* curvature(X) .+ (J / 2π) .* log.(abs.(z .- zc)))
    g = arclength_derivative(X, μ) .* conj.(τ)
    gs = complex.(arclength_derivative(X, real.(g)), arclength_derivative(X, imag.(g)))
    ∇p = C * g .- c .* g .+ gs .* (ds / (2π * im)) .- (J / 2π) ./ (z .- zc)
    -real.(∇p .* (-im .* τ))
end
sol = evolve(θ -> (1 + 0.05cos(3θ) + 0.03sin(5θ) + 0.02cos(8θ)) .* (cos(θ), sin(θ)), hele_shaw, (0, 8.0);
             dt=(X, t) -> 0.5perimeter(X) / length(X), symbol=q -> -0.005q^3, Nmax=1024, saveat=0:0.1:8)
csv("hele_shaw_movie.csv", ("t", "curve", "x", "y"), reduce(vcat, [points_rows(t, 1, X) for (t, X) in zip(sol.t, sol.curves)]))
csv("hele_shaw_area.csv", ("t", "area", "exact", "N"),
    [(t, area(X), area(sol.curves[1]) + 2π * t, length(X)) for (t, X) in zip(sol.t, sol.curves)])

# Hele-Shaw with three bubbles screening each other (examples/hele_shaw_bubbles.jl): one solve for every interface,
# each bubble's pressure an unknown, its source at the deepest point inside it, Δt from the last step's max |U|
function source_point(X)
    x = [complex(p...) for p in curve_points(X)]; n = -im .* cis.(tangent_angles(X)); r = sqrt(area(X) / π)
    inside(q) = abs(sum(angle((x[mod1(j + 1, end)] - q) / (x[j] - q)) for j in eachindex(x))) > π
    cands = filter(inside, [complex(centroid(X)...); x .- 0.25r .* n; x .- 0.5r .* n])
    cands[argmax([minimum(abs.(x .- q)) for q in cands])]
end
function hele_shaw_bubbles(Xs, t; S=0.01, J=fill(2π, length(Xs)))
    nb = length(Xs); e = cumsum(length.(Xs)); blocks = [e[k]-length(Xs[k])+1:e[k] for k in 1:nb]
    z = [complex(x...) for X in Xs for x in curve_points(X)]; τ = reduce(vcat, [cis.(tangent_angles(X)) for X in Xs])
    w = reduce(vcat, [fill(perimeter(X) / length(X), length(X)) for X in Xs]); c = source_point.(Xs); n = length(z)
    C = [i == j ? 0.0im : τ[j] * w[j] / (2π * im * (z[j] - z[i])) for i in 1:n, j in 1:n]
    own = zeros(ComplexF64, n)
    for r in blocks, i in r; own[i] = sum(@view C[i, r]); end
    M = zeros(n + nb, n + nb); M[1:n, 1:n] .= real.(C)
    for i in 1:n; M[i, i] -= real(own[i]); end
    for (k, r) in enumerate(blocks); M[r, n+k] .= -1; M[n+k, r] .= w[r]; end
    f = -S .* reduce(vcat, curvature.(Xs)) .+ sum((J[m] / 2π) .* log.(abs.(z .- c[m])) for m in 1:nb)
    μ = (M \ [f; zeros(nb)])[1:n]
    g = reduce(vcat, [arclength_derivative(Xs[k], μ[r]) .* conj.(τ[r]) for (k, r) in enumerate(blocks)])
    gs = reduce(vcat, [complex.(arclength_derivative(Xs[k], real.(g[r])), arclength_derivative(Xs[k], imag.(g[r])))
                       for (k, r) in enumerate(blocks)])
    ∇p = C * g .- own .* g .+ gs .* w ./ (2π * im) .- sum((J[m] / 2π) ./ (z .- c[m]) for m in 1:nb)
    U = -real.(∇p .* (-im .* τ)); [U[r] for r in blocks]
end
bubble(x0, y0, r, ϕ) = θ -> (x0 + r * (1 + 0.04cos(3θ + ϕ) + 0.03sin(5θ + 2ϕ)) * cos(θ),
                             y0 + r * (1 + 0.04cos(3θ + ϕ) + 0.03sin(5θ + 2ϕ)) * sin(θ))
Xs = [bubble(2.2cos(a), 2.2sin(a), r, a) for (a, r) in zip((π / 2, 7π / 6, 11π / 6), (0.7, 0.5, 0.6))]
umax = Ref(maximum(U -> maximum(abs, U), hele_shaw_bubbles(from_curve.(Xs), 0.0)))
sol = evolve(Xs, hele_shaw_bubbles, (0, 6.0); dt=(Xs, t) -> 0.25minimum(X -> perimeter(X) / length(X), Xs) / umax[],
             symbol=q -> -0.01q^3, Nmax=512, saveat=0:0.1:6,
             callback=integ -> (umax[] = maximum(U -> maximum(abs, U), integ.U);
                                min_distance(integ.states) < 0.05 ? :stop : nothing))
bubble_rows(t, Xs, m) = reduce(vcat, [(x = curve_points(X); k = max(1, length(x) ÷ m);
                                       [(t, id, p[1], p[2]) for p in [x[1:k:end]; [x[1]]]]) for (id, X) in enumerate(Xs)])
csv("bubbles_movie.csv", ("t", "curve", "x", "y"), reduce(vcat, [bubble_rows(t, Xs, 256) for (t, Xs) in zip(sol.t, sol.curves)]))
csv("bubbles_area.csv", ("t", "curve", "area", "exact"),
    [(t, k, area(Xs[k]), area(sol.curves[1][k]) + 2π * t) for (t, Xs) in zip(sol.t, sol.curves) for k in eachindex(Xs)])

# vesicles: Willmore flow with area and length held (examples/vesicles.jl), four reduced areas in a 2 × 2 grid
function willmore(X, A0, L0; τ=0.01)
    κ = curvature(X); f = arclength_derivative(X, arclength_derivative(X, κ)) .+ κ .^ 3 ./ 2
    ∫(g) = curve_integral(X, g)
    p, σ = [perimeter(X) -∫(κ); ∫(κ) -∫(κ .^ 2)] \ [(A0 - area(X)) / τ - ∫(f), (L0 - perimeter(X)) / τ - ∫(κ .* f)]
    f .+ p .- σ .* κ
end
reduced_area(X) = 4π * area(X) / perimeter(X)^2
function ellipse(ν, x0, y0)
    shape(r) = θ -> (x0 + cos(θ) / sqrt(r), y0 + sqrt(r) * sin(θ))
    lo, hi = 1e-3, 1.0
    for _ in 1:50
        r = (lo + hi) / 2; reduced_area(from_curve(shape(r))) < ν ? (lo = r) : (hi = r)
    end
    shape((lo + hi) / 2)
end
Xs = [ellipse(ν, x0, y0) for (ν, x0, y0) in ((0.9, 0, 0), (0.75, 5, 0), (0.6, 0, -2.8), (0.45, 5, -2.8))]
A0 = [area(from_curve(X)) for X in Xs]; L0 = [perimeter(from_curve(X)) for X in Xs]
sol = evolve(Xs, (Xs, t) -> willmore.(Xs, A0, L0), (0, 0.5);
             dt=(Xs, t) -> min(1e-3, 0.0125 / maximum(X -> maximum(abs, curvature(X))^4, Xs)), symbol=q -> -q^4,
             saveat=[0; exp10.(range(-4, log10(0.5), length=60))])
vrows(t, Xs, m=512) = reduce(vcat, [(x = curve_points(X); k = max(1, length(x) ÷ m);
                                      [(t, id, p[1], p[2]) for p in [x[1:k:end]; [x[1]]]]) for (id, X) in enumerate(Xs)])
csv("vesicles_movie.csv", ("t", "curve", "x", "y"), reduce(vcat, [vrows(t, Xs, 128) for (t, Xs) in zip(sol.t, sol.curves)]))
csv("vesicles_energy.csv", ("t", "curve", "nu", "energy"),
    [(t, k, reduced_area(Xs[k]), curve_integral(Xs[k], curvature(Xs[k]) .^ 2) / 2) for (t, Xs) in zip(sol.t, sol.curves) for k in 1:4])

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
