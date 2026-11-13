# Step-size control with ThetaL: surface diffusion of a four-lobed flower, U = κ_ss, with the step chosen from the
# embedded error estimate of ETDRK2. The first steps must be small (the initial layer of the stiff symbol σ = −q⁴:
# the small scales decay within a step and the estimate is first order there), later ones can be large. A step
# whose estimate exceeds `errtol` is rejected before its end evaluation, so a rejection costs one velocity
# evaluation. Drawn: the accepted steps (log₁₀), then the curve.
#   julia --project=examples examples/step_control.jl                   (from the package root)
using ThetaL, DotPlots

X = ξ -> (1 + 0.25cos(4ξ)) .* (cos(ξ), sin(ξ))
κss(X) = arclength_derivative(X, arclength_derivative(X, curvature(X)))
integ = init(CurveProblem([X], (Xs, ev; warm) -> κss.(Xs); symbol=q -> -q^4), ETDRK2())

# a few lines of control: the next step from the estimate, accepted or not (order 2 against order 1: err ~ h²)
function run_controlled!(integ; T, h, tol)
    hs = Float64[]; frames = [[curve_points(integ.states[1])]]; times = [0.0]
    while integ.t < T
        t_before = integ.t
        status = advance!(integ, min(integ.t + h, T); errtol=tol)
        status === :success && push!(hs, integ.t - t_before)
        h *= clamp(0.9 * sqrt(tol / integ.stats.error_estimate), 0.2, 5.0)
        if integ.t >= times[end] + T / 8
            push!(frames, [curve_points(integ.states[1])]); push!(times, integ.t)
            println(progress_line(integ))
        end
    end
    hs, frames, times
end
hs, frames, times = run_controlled!(integ; T=0.02, h=1e-8, tol=1e-4)

dotplot(log10.(hs); height=12)                            # log10 of the accepted steps, one after another
dotplot(frames; t=times, animate=true)
