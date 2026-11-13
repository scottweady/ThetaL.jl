# Mean-curvature flow with ThetaL: a five-lobed star moves with outward normal velocity U = −κ. It turns convex,
# then round, and vanishes at T = A₀/2π. The stiff part of U, the growth rate σ = −q² of a mode of wavenumber
# q = 2πk/L, goes into the exponential step through `symbol`.
#   julia --project=examples examples/mean_curvature_flow.jl            (from the package root)
using ThetaL, DotPlots

X = ξ -> (1 + 0.3cos(5ξ)) .* (cos(ξ), sin(ξ))
T = area(from_curve(X)) / 2π                                         # extinction time

sol = evolve(X, (X, t) -> -curvature(X), (0, 0.95T);
             dt=(X, t) -> 0.0125 / maximum(abs2, curvature(X)),      # Δt ∝ 1/κ²
             symbol=q -> -q^2, saveat=0.95T / 8, verbose=true)
dotplot(sol; animate=true)
