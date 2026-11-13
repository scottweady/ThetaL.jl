# Wind-up in a vortex with ThetaL: a disk off the center of the Lamb–Oseen vortex u = ω(r)(−y, x),
# ω(r) = (1 − e^{−r²})/r², is sheared into a spiral. Each point keeps its radius and turns at ω(r); the curve
# lengthens and sharpens, and N doubles as it needs to. No `symbol`: the motion is not stiff; the step integrates
# the rigid part of the motion (the mean rotation and the drift of the disk) exactly, and the shear explicitly.
#   julia --project=examples examples/vortex_windup.jl                  (from the package root)
using ThetaL, DotPlots

ω(r) = r < 1e-8 ? 1 - r^2 / 2 : -expm1(-r^2) / r^2
X = ξ -> (1.0 + 0.5cos(ξ), 0.5sin(ξ))
U(X, t) = velocity_components(X, [ω(hypot(x...)) .* (-x[2], x[1]) for x in curve_points(X)])[1]

T = 10.0
sol = evolve(X, U, (0, T); dt=(X, t) -> 0.5perimeter(X) / length(X),   # Δt ∝ arclength spacing / |u|, |u| < 1
             saveat=T / 20, verbose=true)
dotplot(sol; animate=true)
