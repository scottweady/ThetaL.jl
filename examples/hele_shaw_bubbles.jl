# Several bubbles in a Hele-Shaw cell with ThetaL: air injected into each of three bubbles at the rate J. Each
# grows fingers, as in hele_shaw.jl, but the bubbles now share the viscous fluid around them: the fluid between
# two bubbles is pushed out by both, its pressure stays high, and the facing sides flatten and hold back while the
# fingers grow outward, into open fluid. The bubbles screen each other.
#
# Outside the bubbles p is harmonic, with u = −∇p and p ~ −(ΣJ/2π) log r at infinity; on bubble k's interface
# p = P_k − Sκ, with P_k the pressure inside it, one unknown constant per bubble, set by its injection rate:
# ∫U ds = J over its interface. U = −∂p/∂n. Several curves move together under one coupled velocity U(Xs, t),
# one vector per curve; the stiff part is the same as for one bubble, σ ≈ −S q³, curve by curve.
#   julia --project=examples examples/hele_shaw_bubbles.jl              (from the package root)
using ThetaL, DotPlots, LinearAlgebra

const S = 0.01                   # surface tension (twice hele_shaw.jl's: fewer, broader fingers, and N ≤ 512)
const J = [2π, 2π, 2π]           # injection rate into each bubble: its area grows at J[k]

# The velocity. As for one bubble (hele_shaw.jl), p = Σ_k −(J_k/2π) log|z − c_k| + Re Φ with c_k inside bubble k,
# and Φ a Cauchy integral of a real density μ, now over every interface. The logs carry each bubble's flux, so Φ
# carries none, and Re Φ → 0 at infinity fixes the constant. At a point z of curve k the singularity is
# subtracted on curve k only (the other curves are smooth there):
#     Φ(z) = (1/2πi)∮_k (μ(ζ) − μ(z))/(ζ − z) dζ + Σ_{j≠k} (1/2πi)∮_j μ(ζ)/(ζ − z) dζ.
# A density constant on one curve gives Φ = 0 outside it, so each curve's mean of μ is free: ∫_k μ ds = 0 takes
# it out, one row per bubble, and the pressures P_k are the matching columns:
#     Re Φ − P_k = −Sκ + Σ_m (J_m/2π) log|z − c_m|  on curve k,        ∫_k μ ds = 0.
#
# The source point c_k may be any point inside bubble k, but the closer it is to the interface the harder the logs
# are to resolve there. The centroid will not do: pushed apart by their neighbors, the bubbles turn lopsided, and a
# fjord can grow onto the centroid. So the source is the deepest of a few candidates: the centroid, and the grid
# points moved inward by a quarter and a half of the bubble's radius, the one inside (winding number) farthest
# from the interface.
function source_point(X)
    x = [complex(p...) for p in curve_points(X)]
    n = -im .* cis.(tangent_angles(X))                # outward normal
    r = sqrt(area(X) / π)
    inside(q) = abs(sum(angle((x[mod1(j + 1, end)] - q) / (x[j] - q)) for j in eachindex(x))) > π
    depth(q) = minimum(abs(xj - q) for xj in x)
    cands = filter(inside, [complex(centroid(X)...); x .- 0.25r .* n; x .- 0.5r .* n])
    cands[argmax(depth.(cands))]
end

function hele_shaw(Xs, t)
    nb = length(Xs)
    blocks = let e = cumsum(length.(Xs)); [e[k]-length(Xs[k])+1:e[k] for k in 1:nb] end
    z = [complex(x...) for X in Xs for x in curve_points(X)]
    τ = reduce(vcat, [cis.(tangent_angles(X)) for X in Xs])
    w = reduce(vcat, [fill(perimeter(X) / length(X), length(X)) for X in Xs])   # trapezoidal weights, ds
    c = source_point.(Xs)
    n = length(z)

    # the trapezoidal Cauchy kernel over all the curves, and the row sums over each point's own curve
    C = zeros(ComplexF64, n, n)
    for j in 1:n
        a = τ[j] * w[j] / (2π * im)
        for i in 1:n
            i == j || (C[i, j] = a / (z[j] - z[i]))
        end
    end
    own = zeros(ComplexF64, n)
    for r in blocks, i in r
        own[i] = sum(@view C[i, r])
    end

    # the bordered system for (μ, P)
    M = zeros(n + nb, n + nb)
    M[1:n, 1:n] .= real.(C)
    for i in 1:n
        M[i, i] -= real(own[i])
    end
    for (k, r) in enumerate(blocks)
        M[r, n+k] .= -1                              # − P_k on curve k
        M[n+k, r] .= w[r]                            # ∫_k μ ds = 0
    end
    f = -S .* reduce(vcat, curvature.(Xs)) .+ sum((J[m] / 2π) .* log.(abs.(z .- c[m])) for m in 1:nb)
    μ = (M \ [f; zeros(nb)])[1:n]

    # Φ' from outside, the same way: g = μ_s/τ, own singularity subtracted, its diagonal term g_s ds/(2πi)
    g = reduce(vcat, [arclength_derivative(Xs[k], μ[r]) .* conj.(τ[r]) for (k, r) in enumerate(blocks)])
    gs = reduce(vcat, [complex.(arclength_derivative(Xs[k], real.(g[r])), arclength_derivative(Xs[k], imag.(g[r])))
                       for (k, r) in enumerate(blocks)])
    dΦ = C * g .- own .* g .+ gs .* w ./ (2π * im)

    ∇p = dΦ .- sum((J[m] / 2π) ./ (z .- c[m]) for m in 1:nb)     # p_x − i p_y
    U = -real.(∇p .* (-im .* τ))
    [U[r] for r in blocks]
end

# three bubbles of different sizes, each a little off round
bubble(x0, y0, r, ϕ) = ξ -> (x0 + r * (1 + 0.04cos(3ξ + ϕ) + 0.03sin(5ξ + 2ϕ)) * cos(ξ),
                             y0 + r * (1 + 0.04cos(3ξ + ϕ) + 0.03sin(5ξ + 2ϕ)) * sin(ξ))
Xs = [bubble(2.2cos(a), 2.2sin(a), r, a) for (a, r) in zip((π / 2, 7π / 6, 11π / 6), (0.7, 0.5, 0.6))]

# Δt: the stiff part is in the symbol, so the step only has to follow the points, Δt ∝ spacing/max|U|. The
# fastest points are the finger tips, so max|U| is read off the velocities of the last step (integ.U), kept by
# the callback. The callback also stops the run when two interfaces (or two fingers of one) come within `gap`:
# below that the coupling between them is stiff in a way no symbol expresses (see the README's Limits), and the
# bubbles would soon merge, which ThetaL does not do.
T, gap, cfl = 6.0, 0.05, 0.25
umax = Ref(maximum(U -> maximum(abs, U), hele_shaw(from_curve.(Xs), 0.0)))
function watch(integ)
    umax[] = maximum(U -> maximum(abs, U), integ.U)
    min_distance(integ.states) < gap ? :stop : nothing
end
sol = evolve(Xs, hele_shaw, (0, T);
             dt=(Xs, t) -> cfl * minimum(X -> perimeter(X) / length(X), Xs) / umax[], symbol=q -> -S * q^3,
             callback=watch, tol=1e-6, Nmax=512, saveat=T / 12, verbose=true)

# each bubble's area grows at its own J: a check on the velocity and on the pressures P_k
A0 = area.(sol.curves[1])
for (t, Xs) in zip(sol.t, sol.curves)
    println("t = ", rpad(round(t; digits=2), 5), "  N = ", rpad(join(length.(Xs), "/"), 12),
            "  area − (A₀ + Jt) = ", join((round(x; sigdigits=2) for x in area.(Xs) .- A0 .- J .* t), ", "))
end
dotplot(sol; animate=true)
