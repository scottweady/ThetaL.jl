# Coarsening with ThetaL: twelve particles coupled through a mean field, U_i = κ̄ − κ_i with κ̄ = Σ∫κ ds / ΣL. Small
# particles shrink and vanish, large ones grow, and the total area is kept. Shows a coupled velocity (every
# particle's U depends on all of them), the integrator driven by hand (`init`, `advance!`), particles removed as
# they vanish (`remove_curve!`, which makes the next step evaluate the velocities afresh), and a step rejected when
# a particle would vanish inside it.
#   julia --project=examples examples/coarsening.jl                     (from the package root)
using ThetaL, DotPlots, Random

Random.seed!(3)
particles = [((3.0i, 3.0j), 0.3 + 0.7rand()) for i in 1:4 for j in 1:3]
Xs = [ξ -> (xc[1] + r * cos(ξ), xc[2] + r * sin(ξ)) for (xc, r) in particles]

# the mean field; remove_area_rate makes Σ∫U ds = 0 exact on the grid (one constant over all particles)
function velocity(Xs, ev; warm)
    κbar = 2π * length(Xs) / sum(perimeter, Xs)
    remove_area_rate(Xs, [κbar .- curvature(X) for X in Xs]; per_curve=false)
end

integ = init(CurveProblem(Xs, velocity; symbol=q -> -q^2), ETDRK2())
radius(X) = perimeter(X) / 2π
frames = [curve_points.(integ.states)]; times = [0.0]

while integ.t < 3 && length(integ.states) > 1
    h = 0.05 * minimum(radius, integ.states)^2              # the smallest particle sets the time scale
    if advance!(integ, integ.t + h) === :rejected           # a particle would reach zero length inside the step
        remove_curve!(integ, argmin(perimeter.(integ.states)))
        continue
    end
    for i in reverse(eachindex(integ.states))               # particles below a size are taken out
        radius(integ.states[i]) < 0.05 && remove_curve!(integ, i)
    end
    if integ.t >= times[end] + 0.5
        push!(frames, curve_points.(integ.states)); push!(times, integ.t)
        println(progress_line(integ))
    end
end
push!(frames, curve_points.(integ.states)); push!(times, integ.t)
dotplot(frames; t=times, animate=true)
