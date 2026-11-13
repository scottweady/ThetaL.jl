# Helpers shared by the test files: test curves and the distance from points to a curve.

# a closed curve given by its tangent angle: θ = α + φ with φ of even harmonics only (closed by symmetry)
smooth_state(N) = ThetaLState([0.3sin(2α) + 0.1cos(4α) for α in 2π .* (0:N-1) ./ N], 2π, Vec2(0.2, -1.0))

# the distance from the points xs to the curve z(α) (dense search, then Newton on (z(α) − x)·z′(α) = 0)
function curve_distance(z, xs)
    αs = range(0, 2π, length=4097)[1:end-1]; zs = z.(αs); h = 1e-4
    maximum(xs) do x
        a = αs[argmin([norm(y - x) for y in zs])]
        for _ in 1:8
            d1 = (z(a + h) - z(a - h)) / 2h; d2 = (z(a + h) - 2z(a) + z(a - h)) / h^2
            a -= dot(z(a) - x, d1) / (dot(d1, d1) + dot(z(a) - x, d2))
        end
        norm(z(a) - x)
    end
end

# `evolve`'s U, or a problem's velocity, for a velocity field u(x): its normal part at the points of the curve(s)
field(u) = (X, t) -> velocity_components(X, u.(curve_points(X)))[1]
field_velocity(u) = (Xs, ev; warm) -> [velocity_components(X, u.(curve_points(X)))[1] for X in Xs]
