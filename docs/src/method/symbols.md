# Symbols, filter and resolution

## The stiff symbol

The small-scale decomposition of Hou, Lowengrub & Shelley: at high wavenumbers a velocity built from the curve's
own geometry relaxes the tangent angle at a rate that grows with the wavenumber,
``\theta_t \approx \sigma_k\hat\phi_k`` with ``\sigma_k < 0``. The caller gives that growth rate as a function of
the wavenumber along the curve, ``q = 2\pi|k|/L``: `symbol = q -> σ`, or `symbol = (q, X) -> σ` when it depends on
the curve `X`. It is the dispersion relation of the velocity linearized about a straight interface. The step
evaluates it once per curve at the start of each step and holds the rates for the step. Whatever part of the
velocity the symbol does not capture stays in the explicit remainder, so the symbol need not be exact: the step
converges at its order with a wrong symbol too ([Method overview](@ref), panel (d)), but steps that are large
compared with ``1/|\sigma|`` need the stiff part captured.

| Velocity | Symbol |
|:--|:--|
| curve shortening, ``U = -\kappa`` | `q -> -q^2` |
| Hele-Shaw with surface tension ``S`` (sheet strength ``S\kappa_s - 2A_\mu\,u\cdot\tau``, as in HLS 1994) | `q -> -S / 2 * q^3` |
| surface diffusion, ``U = \kappa_{ss}`` | `q -> -q^4` |
| anisotropic, ``U = -(\gamma + \gamma'')(\vartheta)\,\kappa`` | `(q, X) -> -maximum(γ + γ″ on X) * q^2` |
| a numerical symbol (a table, a quadrature per mode) | any function of ``q`` returning one rate |

A symbol that depends on the curve must bound the stiffness it stands for (anisotropy: the largest value of
``\gamma + \gamma''`` on the curve, `examples/anisotropy.jl`). A positive ``\sigma`` makes the exponential amplify
that mode; ThetaL warns. A rate that is not a finite real number is an error.

## The spectral filter

Every scheme carries a `SpectralFilter(; order = 36, cutoff = 0)`: mode ``k`` is multiplied by
``\exp(-36\,(|k|/k_c)^{\mathrm{order}})`` and modes ``|k| \ge k_c`` are removed, on ``\phi`` after every stage and on
the deformation part of ``U`` before it is differentiated. `cutoff = 0` means ``k_c = N/2``; an integer sets it;
`:auto` sets it once, after the first velocity evaluation, from `noise_cutoff`: a computed velocity carries an
error about flat in ``k``, which the step damps at the rate of the symbol, so the noise settles at
``\hat\phi_k \approx (2\pi/L) k \varepsilon/|\sigma_k|``; the cutoff is placed where the curve's spectrum meets that
floor. Set it through the scheme: `ETDRK2(filter_cutoff = :auto)`.

## Resolution

The tail of a curve's spectrum (`spectral_tail`: the largest ``|\hat\phi_k|/N`` around ``k = N/2``) measures how well
N resolves it. `from_curve` picks the smallest power of two from 64 whose tail is below `tol/100`. `evolve` and
`solve` double N whenever the tail exceeds `tol` (up to `Nmax`), and, with `coarsen = true`, halve it once the tail
on half the points has stayed below `tol/100` for ten steps; the gap between the two thresholds and the wait keep
them from undoing each other. One step over `tol` is enough to double N, and without `coarsen` it stays doubled, so
a computed velocity whose error spikes now and then can push N up for good: `ETDRK2(filter_cutoff = :auto)` damps
that noise (above), and `coarsen = true` brings N back. `verbose = true` prints each change of N with its time and
cause. The default `Nmax = 16384` is far beyond what an O(N²) velocity solve can afford; set it. With the integrator, `double_grid!` and `halve_grid!` change N between steps (the
velocities are interpolated or truncated with the curve), on whatever rule the caller keeps.
