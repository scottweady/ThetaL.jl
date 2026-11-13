# Anisotropic curve shortening with ThetaL: U = −(γ + γ″)(ϑ) κ with the surface energy γ(ϑ) = 1 + ε cos 4ϑ, ϑ the
# angle of the normal. A circle shrinks toward the four-fold Wulff shape (drawn scaled to area π). The stiffness
# (γ + γ″)(ϑ) varies along the curve, so the symbol σ(q, X) reads it from the curve X: frozen per curve and per
# step at its largest value (a bound, as a frozen symbol must be).
#   julia --project=examples examples/anisotropy.jl                     (from the package root)
using ThetaL, DotPlots

ε = 0.05                                                    # γ + γ″ = 1 − 15ε cos 4ϑ > 0 for ε < 1/15
stiffness(ϑ) = 1 - 15ε * cos(4ϑ)
normal_angles(X) = tangent_angles(X) .- π / 2

velocity(Xs, ev; warm) = [-stiffness.(normal_angles(X)) .* curvature(X) for X in Xs]
symbol(q, X) = -maximum(stiffness, normal_angles(X)) * q^2

X = ξ -> (cos(ξ), sin(ξ))
T = 0.45                                                     # the circle of area π would vanish at t = 0.5
sol = solve(CurveProblem([X], velocity; symbol), ETDRK2(); tspan=(0, T),
            dt=(Xs, t) -> 0.01 / maximum(abs2, curvature(Xs[1])), saveat=T / 9, verbose=true)
dotplot([[curve_points(Xs[1]) ./ sqrt(area(Xs[1]) / π)] for Xs in sol.curves]; t=sol.t, animate=true)
