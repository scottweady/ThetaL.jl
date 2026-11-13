# Vesicles with ThetaL: closed membranes that resist bending (Helfrich), whose enclosed area A and length L are
# both fixed (the fluid inside cannot leave, the membrane cannot stretch). Each relaxes by Willmore flow, the
# gradient flow of the bending energy E = ½∮κ² ds, to the shape set by its reduced area ν = 4πA/L² (1 for a
# circle): an ellipse when ν is near 1, then a biconcave shape, the two-dimensional red blood cell, and at low ν a
# deeper and deeper dimple.
#
# With U the outward normal velocity, κ_t = −U_ss − κ²U and (ds)_t = κU ds, so dE/dt = −∮(κ_ss + ½κ³) U ds: the flow
# is U = κ_ss + ½κ³, and with the constraints
#     U = κ_ss + ½κ³ + p − σκ,
# a pressure p and a tension σ chosen at every evaluation to make dA/dt = ∮U ds and dL/dt = ∮κU ds vanish (a 2×2
# system per curve). Those rates are kept to zero, but A and L would still drift by the error of the step, most
# in the first steps, while the small scales of the starting shape decay faster than a step (the initial layer of
# the README's Limits); so the rates are set instead to (A₀ − A)/τ and (L₀ − L)/τ, which pulls A and L back to
# their targets. The stiff part is κ_ss, so the symbol is −q⁴, as for surface diffusion. The curves do not
# interact: four vesicles in one call, U(Xs, t) one vector per curve.
#   julia --project=examples examples/vesicles.jl                       (from the package root)
using ThetaL, DotPlots

function willmore(X, A0, L0; τ=0.01)
    κ = curvature(X)
    f = arclength_derivative(X, arclength_derivative(X, κ)) .+ κ .^ 3 ./ 2
    ∫(g) = curve_integral(X, g)
    # ∫(f + p − σκ) ds = (A₀ − A)/τ and ∫κ(f + p − σκ) ds = (L₀ − L)/τ for (p, σ)
    p, σ = [perimeter(X) -∫(κ); ∫(κ) -∫(κ .^ 2)] \
           [(A0 - area(X)) / τ - ∫(f), (L0 - perimeter(X)) / τ - ∫(κ .* f)]
    f .+ p .- σ .* κ
end
bending(X) = curve_integral(X, curvature(X) .^ 2) / 2
reduced_area(X) = 4π * area(X) / perimeter(X)^2

# ellipses of area π with the reduced areas ν, stacked (their aspect ratio by bisection on ν)
function ellipse(ν, y0)
    shape(r) = ξ -> (cos(ξ) / sqrt(r), y0 + sqrt(r) * sin(ξ))         # semi-axes 1/√r, √r: area π
    lo, hi = 1e-3, 1.0
    for _ in 1:50
        r = (lo + hi) / 2
        reduced_area(from_curve(shape(r))) < ν ? (lo = r) : (hi = r)
    end
    shape((lo + hi) / 2)
end
νs = [0.9, 0.75, 0.6, 0.45]
Xs = [ellipse(ν, -2.5(k - 1)) for (k, ν) in enumerate(νs)]
A0 = [area(from_curve(X)) for X in Xs]
L0 = [perimeter(from_curve(X)) for X in Xs]

# Δt ∝ 1/max κ⁴ (the explicit terms, κ³ and the constraints, are as stiff there as κ_ss): small steps while the
# tips of the ellipses relax, then up to 10⁻³. Larger first steps keep A and L (the feedback sees to that) but
# leave an error that breaks the symmetry of the start: the most deflated vesicle drifts sideways.
T = 0.5
sol = evolve(Xs, (Xs, t) -> willmore.(Xs, A0, L0), (0, T);
             dt=(Xs, t) -> min(1e-3, 0.0125 / maximum(X -> maximum(abs, curvature(X))^4, Xs)),
             symbol=q -> -q^4, saveat=[0, 0.003, 0.01, 0.03, 0.1, T], verbose=true)

# the energy falls, and ν holds
for (t, Xs) in zip(sol.t, sol.curves)
    println("t = ", rpad(t, 5), "  E = ", join((rpad(round(bending(X); digits=4), 8) for X in Xs), " "),
            "  ν − ν₀ = ", join((round(reduced_area(X) - ν; sigdigits=2) for (X, ν) in zip(Xs, νs)), ", "))
end
dotplot(sol; animate=true)
