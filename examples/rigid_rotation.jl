# Rigid rotation with ThetaL: a star off the origin turns once about it in the flow u = ω(−y, x), moved by its
# normal velocity U = u·n alone. The step takes the rigid motion out of U and integrates it exactly, so the star
# keeps its shape at any Δt, here ten steps a turn.
#   julia --project=examples examples/rigid_rotation.jl                 (from the package root)
using ThetaL, DotPlots

ω, xc = 2π, (1.0, 0.0)                                                # one turn in t = 1
X = ξ -> xc .+ 0.5 * (1 + 0.3cos(5ξ)) .* (cos(ξ), sin(ξ))
U(X, t) = velocity_components(X, [ω .* (-x[2], x[1]) for x in curve_points(X)])[1]     # U = u·n

sol = evolve(X, U, (0, 1); dt=0.1, saveat=0.1, verbose=true)
dotplot(sol; animate=true)
