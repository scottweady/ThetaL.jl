# Hele-Shaw with ThetaL: air injected at rate J into a bubble surrounded by viscous fluid. The bubble grows, and
# where it grows fast enough surface tension S cannot hold it round: it splits into fingers (Saffman–Taylor).
#
# Outside the bubble the pressure p is harmonic, with Darcy's law u = −∇p and p ~ −(J/2π) log r at infinity; on
# the interface p = −Sκ (Laplace's law, up to the bubble pressure, a constant). The interface moves with the
# fluid, U = −∂p/∂n. ThetaL knows none of this: the velocity below solves the exterior Dirichlet problem with a
# boundary integral, and hands back U at the grid points.
#
# The stiff part: on a circle of radius R a mode k grows at σ = (k − 1)J/(2πR²) − S k(k² − 1)/R³, so at small
# scales σ ≈ −S q³ with q = k/R the wavenumber along the curve. That goes into the exponential step as the
# symbol. (HLS 1994 write the same problem with Darcy's law scaled by 1/2, hence their −S/2 q³.)
#   julia --project=examples examples/hele_shaw.jl                      (from the package root)
using ThetaL, DotPlots, LinearAlgebra

const S = 0.005     # surface tension
const J = 2π        # injection rate: the area grows at J

# The velocity. In complex form z = x + iy, write p = −(J/2π) log|z − z_c| + Re Φ, Φ analytic outside the bubble
# (z_c, the centroid, is any point inside), and Φ(z) = (1/2πi)∮ μ(ζ)/(ζ − z) dζ + ∫μ ds, a Cauchy integral of a real
# density μ plus a constant. Its limit from outside, after taking out the singularity, is
#     Φ(z) = (1/2πi)∮ (μ(ζ) − μ(z))/(ζ − z) dζ + ∫μ ds,
# smooth, so the trapezoidal rule on ThetaL's equal-arclength grid is spectrally accurate; on the diagonal the
# integrand is μ_s τ, imaginary for real μ, so its real part drops. Re Φ = −Sκ + (J/2π) log|z − z_c| on the
# boundary is a dense N × N system for μ (the ∫μ ds term fixes the constant that the Cauchy integral misses).
# Then Φ' = (1/2πi)∮ g/(ζ − z) dζ with g = μ_s/τ (by parts), limited from outside the same way, gives
# p_x − i p_y = Φ' − (J/2π)/(z − z_c), and U = −∂p/∂n = −Re[(p_x − i p_y) n], n = −iτ outward.
function hele_shaw(X, t)
    N, ds = length(X), perimeter(X) / length(X)
    z = [complex(x...) for x in curve_points(X)]
    τ = cis.(tangent_angles(X))                       # unit tangent, counterclockwise
    zc = complex(centroid(X)...)

    # the trapezoidal Cauchy kernel C[i, j] = τ_j ds / (2πi (z_j − z_i)), zero on the diagonal, and its row sums
    # (the subtracted singularity: Σ_j C[i, j] (f_j − f_i) = (C f)_i − c_i f_i)
    C = zeros(ComplexF64, N, N)
    for j in 1:N
        w = τ[j] * ds / (2π * im)
        for i in 1:N
            i == j || (C[i, j] = w / (z[j] - z[i]))
        end
    end
    c = vec(sum(C; dims=2))

    # Re Φ = Re C μ − Re c μ + ∫μ ds
    A = real.(C)
    A[diagind(A)] .-= real.(c)
    A .+= ds
    μ = A \ (-S .* curvature(X) .+ (J / 2π) .* log.(abs.(z .- zc)))

    # Φ' from outside: (1/2πi)∮ (g(ζ) − g(z))/(ζ − z) dζ, the diagonal term g_s ds/(2πi)
    g = arclength_derivative(X, μ) .* conj.(τ)
    gs = complex.(arclength_derivative(X, real.(g)), arclength_derivative(X, imag.(g)))
    dΦ = C * g .- c .* g .+ gs .* (ds / (2π * im))

    ∇p = dΦ .- (J / 2π) ./ (z .- zc)                  # p_x − i p_y
    -real.(∇p .* (-im .* τ))
end

X0 = ξ -> (1 + 0.05cos(3ξ) + 0.03sin(5ξ) + 0.02cos(8ξ)) .* (cos(ξ), sin(ξ))
T = 8.0
# Nmax: the dense solve is O(N³). N doubles as the fingers sharpen (256 → 512 → 1024 here). Leave the filter at its
# default: filter_cutoff = :auto would also strip the fine scales the fingers really have, hiding them from the
# refinement test, and the area check below would show it.
# Δt: the surface tension is in the symbol, so the step only has to follow the points as they move along the
# curve, Δt ≲ spacing/|u| with |u| ≲ 1 here; it shrinks as N doubles.
sol = evolve(X0, hele_shaw, (0, T); dt=(X, t) -> 0.5perimeter(X) / length(X), symbol=q -> -S * q^3,
             tol=1e-6, Nmax=1024, saveat=T / 8, verbose=true)

# the area grows at exactly J: a check on the velocity, since ∫U ds = J is not imposed anywhere
for (t, X) in zip(sol.t, sol.curves)
    println("t = ", rpad(round(t; digits=2), 5), "  N = ", rpad(length(X), 5),
            "  area − (A₀ + Jt) = ", area(X) - area(sol.curves[1]) - J * t)
end
dotplot(sol; animate=true)
