# Step-size control with ThetaL: surface diffusion of a four-lobed flower, U = κ_ss, with the step chosen from the
# embedded error estimate of ETDRK2. The first steps must be small (the initial layer of the stiff symbol σ = −q⁴:
# the small scales decay within a step and the estimate is first order there), later ones can be large. A step
# whose estimate exceeds `errtol` is rejected before its end evaluation, so a rejection costs one velocity
# evaluation. In one call, `evolve(…; adaptive = true, errtol)`; by hand, the same few lines of control around
# `advance!`. Drawn: the accepted steps (log₁₀), then the curve.
#   julia --project=examples examples/step_control.jl                   (from the package root)
using ThetaL, DotPlots

X = ξ -> (1 + 0.25cos(4ξ)) .* (cos(ξ), sin(ξ))
κss(X) = arclength_derivative(X, arclength_derivative(X, curvature(X)))

# in one call: the driver picks the steps (and lands on the saved times); the callback records them
hs = Float64[]
sol = evolve(X, (X, t) -> κss(X), (0, 0.02); symbol=q -> -q^4, adaptive=true, errtol=1e-4, saveat=0.02 / 8,
             callback=integ -> (push!(hs, integ.stats.last_step); nothing), verbose=true)
println("evolve: ", sol.steps, " steps, ", sol.rejections, " rejected, from Δt = ", hs[1], " to ", maximum(hs))

# by hand: the next step from the estimate, accepted or not (order 2 against order 1: err ~ h²)
integ = init(CurveProblem([X], (Xs, ev; warm) -> κss.(Xs); symbol=q -> -q^4), ETDRK2())
function run_controlled!(integ; T, h, tol)
    hs = Float64[]
    while integ.t < T
        t_before = integ.t
        status = advance!(integ, min(integ.t + h, T); errtol=tol)
        status === :success && push!(hs, integ.t - t_before)
        h *= clamp(0.9 * sqrt(tol / integ.stats.error_estimate), 0.2, 5.0)
    end
    hs
end
hs_hand = run_controlled!(integ; T=0.02, h=1e-8, tol=1e-4)
println("by hand: ", progress_line(integ))

dotplot(log10.(hs); height=12)                            # log10 of the accepted steps, one after another
dotplot(sol; animate=true)
