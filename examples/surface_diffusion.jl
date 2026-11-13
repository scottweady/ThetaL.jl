# Surface diffusion with ThetaL: a four-lobed flower moves with outward normal velocity U = κ_ss and relaxes to the
# circle of the same area; its lobes decay at the rate q⁴ ≈ 256, so the run is short. The stiff part of U, the
# growth rate σ = −q⁴ of a mode of wavenumber q = 2πk/L, goes into the exponential step through `symbol`, so the
# step is not tied to Δt ∝ (L/N)⁴.
#   julia --project=examples examples/surface_diffusion.jl              (from the package root)
using ThetaL, DotPlots

X = ξ -> (1 + 0.25cos(4ξ)) .* (cos(ξ), sin(ξ))
κss(X) = arclength_derivative(X, arclength_derivative(X, curvature(X)))

sol = evolve(X, (X, t) -> κss(X), (0, 0.02); dt=2.5e-5, symbol=q -> -q^4, saveat=0.02 / 10,
             verbose=true)
dotplot(sol; animate=true)
