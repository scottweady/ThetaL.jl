# Geometry and conservation utilities for callers: distances between and within curves (the close-approach monitor),
# integrals along a curve, the parallel curve, and the projection of a normal velocity that removes area change.

"""
    min_distance(states; within = true) -> Float64

The smallest distance between grid points of different curves, and, with `within`, between points of one curve
that are more than a quarter of its length apart along it (a neck). The small-scale decomposition fails where this
becomes comparable to the grid spacing; use it to stop, refine or shorten the step. O(N²) per pair of curves.
"""
function min_distance(states::AbstractVector{ThetaLState}; within::Bool=true)
    pts = [curve_points(X) for X in states]
    d = Inf
    for a in eachindex(pts), b in a+1:length(pts), x in pts[a], y in pts[b]
        d = min(d, hypot(x[1] - y[1], x[2] - y[2]))
    end
    if within
        for p in pts
            N = length(p)
            for i in 1:N, j in i+1:N
                min(j - i, N - (j - i)) > N ÷ 4 || continue
                d = min(d, hypot(p[i][1] - p[j][1], p[i][2] - p[j][2]))
            end
        end
    end
    d
end

"""
    curve_integral(X, f) -> ∫ f ds

The integral along the curve `X` of `f` given at its grid points (equal arclength: the trapezoidal rule, spectrally
accurate for smooth periodic `f`).
"""
function curve_integral(X::ThetaLState, f::AbstractVector)
    length(f) == length(X.φ) || throw(DimensionMismatch("$(length(f)) values for $(length(X.φ)) points"))
    sum(f) * X.L / length(f)
end

"""
    offset_curve(X, δ) -> ThetaLState

The parallel curve at distance `δ` along the outward normal (δ < 0: inward), on the same number of points: the
points x(α) + δ n(α), put back at equal arclength (its length is L + 2πδ, its area A + δL + πδ²). Valid while
δ > −1/max κ (the parallel curve of a convex part stays smooth only so far inward).
"""
function offset_curve(X::ThetaLState, δ::Real)
    X.L + 2π * δ > 0 || throw(ArgumentError("offset $δ leaves no curve (length $(X.L + 2π * δ))"))
    z = curve_function(X); φ = _trig_interpolant(X.φ)
    parallel(α) = (θ = α + φ(α); z(α) + δ * Vec2(sin(θ), -cos(θ)))
    from_curve(parallel, length(X.φ))
end

"""
    remove_area_rate(states, U; per_curve = true) -> U′

The normal velocities with the part that changes area taken out: ∫ U′ ds = 0 on every curve (`per_curve = true`), or
one constant subtracted from all curves so that only the total area is kept, Σᵢ ∫ U′ᵢ ds = 0 (`per_curve = false`,
e.g. Ostwald ripening, where area moves between particles). Area then drifts only at the time error of the step.
"""
function remove_area_rate(states::AbstractVector{ThetaLState}, U::AbstractVector; per_curve::Bool=true)
    means = [sum(u) / length(u) for u in U]                  # ∫ U ds / L on each curve (equal arclength)
    if per_curve
        return [Vector{Float64}(u) .- m for (u, m) in zip(U, means)]
    end
    shift = sum(X.L * m for (X, m) in zip(states, means)) / sum(X.L for X in states)
    [Vector{Float64}(u) .- shift for u in U]
end
