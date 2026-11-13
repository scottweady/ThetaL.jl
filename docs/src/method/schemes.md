# Schemes

## One step of all curves

At small scales the tangent angle of each curve relaxes as ``\theta_t \approx \sigma_k \hat\phi_k``, with the
growth rate ``\sigma_k = \sigma(2\pi|k|/L) < 0`` of mode ``k`` given by the symbol ([Symbols, filter and
resolution](@ref)). Each curve's equation is written as

```math
\hat\phi_t = \sigma_k \hat\phi_k + \hat R_k, \qquad R = \theta_t - \sigma \hat\phi,
```

with the remainder ``R`` treated explicitly. A one-step exponential scheme takes, per curve and mode, with
``z = \sigma_k h``:

```math
\hat\phi^{(i)} = e^{c_i z} \hat\phi + h \sum_{j < i} A_{ij}(z)\, \hat R^{(j)},
```

```math
\hat\phi_{n+1} = e^{z} \hat\phi + h \sum_i b_i(z)\, \hat R^{(i)}.
```

The length uses the same weights at ``z = 0``, and ``x_0`` the rigid motions of the stages averaged with them (exact
when the motion is constant over the step). Each stage after the first evaluates the velocity once, for all curves
together: the curves are coupled only through it. The velocity at the end of a step is the next step's first stage.

## The schemes

| Scheme | ``c`` | ``A_{21}(z)`` | ``b(z)`` | Order | Evaluations per step | Estimate |
|:--|:--|:--|:--|:--|:--|:--|
| `ETDRK2()` (default) | ``0, \tfrac12`` | ``\tfrac12 \varphi_1(z/2)`` | ``\varphi_1 - 2\varphi_2,\ 2\varphi_2`` | 2 | 2 | embedded `ExpEuler` |
| `ExpEuler()` | ``0`` | | ``\varphi_1`` | 1 | 1 | none |

with ``\varphi_1(z) = (e^z - 1)/z`` and ``\varphi_2(z) = (e^z - 1 - z)/z^2`` (by their series near 0). ETDRK2 is the
exponential Runge–Kutta method of Hochbruck & Ostermann (2005) with ``c_2 = 1/2``. With no symbol (``z = 0``) they
are the explicit midpoint rule and Euler's method.

Why exponential weights rather than an integrating factor: for a state at rest under a symbol (``\theta_t = 0``,
so ``R = -\sigma\hat\phi``) the exponential step returns ``e^{z}\hat\phi - h\varphi_1(z)\sigma\hat\phi = \hat\phi``,
whatever ``h``. An integrating-factor (Lawson) step does not, so a shape that only turns would drift. Lawson's
midpoint rule remains in the code as `ThetaL.Lawson2`, a reference for the tests.

## The error estimate

ETDRK2's end state and its embedded exponential Euler solution differ by ``2h\varphi_2(z)(\hat R^{(2)} - \hat
R^{(1)})``, known after the midpoint stage, before the end evaluation. `advance!(integ, t; errtol)` rejects the step
there when the estimate exceeds `errtol`, so a rejected step costs one velocity evaluation. The estimate of the
last step is `integ.stats.error_estimate` (relative, the largest over the curves of ``\phi``, ``L`` and ``x_0``).
It is ``O(h^2)`` once the small scales have relaxed and ``O(h)`` in an initial layer ([Limits](@ref)); `evolve` and
`solve` control the step on it with `adaptive = true`: the next step is ``h \cdot \mathrm{clamp}(0.9\,
(\mathrm{errtol}/\mathrm{err})^{1/2}, 0.2, 5)`` (the estimate is that of the order-1 step, ``O(h^2)``), a step over
`errtol` is retried with the step so shrunk. `examples/step_control.jl` runs it, and the same few lines by hand.

## Adding a scheme

A scheme is a type with one field, the spectral filter, and three methods (`src/stepping/schemes.jl`):

```julia
struct MyScheme <: ThetaL.OneStepScheme
    filter::SpectralFilter
end
MyScheme(; filter_order = 36, filter_cutoff = 0) =
    MyScheme(SpectralFilter(order = filter_order, cutoff = filter_cutoff))
ThetaL.stage_nodes(::MyScheme) = (0.0, 0.5)               # c₁ = 0, c₂, …
ThetaL.scheme_weights(::MyScheme, z::Float64) = (A, b, b̂)  # tuples: A[i][j] (j < i), b[i], b̂[i]
ThetaL.has_error_estimate(::MyScheme) = true               # if b̂ differs from b
```

The generic step (`src/stepping/step.jl`) does the rest: the stages, the velocity evaluations, the rigid motion,
the filter, the estimate, the rejection. Multistep schemes are not supported.
