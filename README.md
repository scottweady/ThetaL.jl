# ThetaL.jl

Boundary evolution of closed plane curves in θ–L variables (Hou, Lowengrub & Shelley 1994). N curves move with
outward normal velocities that may depend on all of them,

    (U₁, …, U_N) = F(X₁, …, X_N; t),

which the caller computes; ThetaL keeps the curves at equal arclength, treats the stiff small scales exactly with
an exponential integrator (from a symbol the caller gives), and integrates the rigid part of each curve's motion
exactly, so a curve that translates and turns without deforming does so exactly at any time step. ThetaL knows
nothing of the physics behind F: in the examples it is curve shortening, surface diffusion, a flow field, or a
mean field coupling many particles.

A curve is held as its tangent angle θ(α) and length L, α ∈ [0, 2π) an equal-arclength parameter, sampled at N
equispaced α (`ThetaLState`: φ = θ − α, L, and the point x₀ at α = 0). Conventions (outward normal
n = (sin θ, −cos θ), κ = θ_s > 0 on convex parts) and the equations are in the module docstring, `?ThetaL`.

Documentation: https://scottweady.github.io/ThetaL.jl/ (getting started, the method, examples, the API reference
and a code map). To build it locally: `julia --project=docs -e 'using Pkg; Pkg.instantiate()'` once, then
`julia --project=docs docs/make.jl`, and open `docs/build/index.html`.

## One call

```julia
using ThetaL, DotPlots

# curve shortening of a five-lobed star: U = −κ, whose stiff part is the growth rate σ = −q² at wavenumber q
X = ξ -> (1 + 0.3cos(5ξ)) .* (cos(ξ), sin(ξ))
sol = evolve(X, (X, t) -> -curvature(X), (0, 0.1); dt = 1e-3, symbol = q -> -q^2, saveat = 0.02)

sol.t                                # 0, 0.02, …, 0.1 (landed on exactly)
area(sol(0.1))                       # ≈ A₀ − 2π·0.1
save_curves("star.csv", sol)         # t, curve id, x, y for every saved curve
dotplot(sol)                         # the saved curves in the terminal, colored by time (DotPlots)
```

`evolve(curve, U, tspan; dt, symbol, saveat, tol, Nmax, scheme, rigid, callback, coarsen, verbose, adaptive, errtol,
dtmin, dtmax)`: `curve` is a parametrization `ξ -> (x, y)`, a `ThetaLState`, or a vector of them; `U(X, t)` (or
`U(Xs, t)` for several) returns the normal velocity at the grid points; `dt` is a number or `dt(X, t)`, or with
`adaptive = true` the step is chosen from ETDRK2's embedded error estimate, kept below `errtol` (default 1e-5;
`dt` then the first step, or a bound if a function; `dtmin`, `dtmax`); N is chosen to `tol`, doubled when needed
and, with `coarsen = true`, halved when the curve has become smooth; `callback(integ)` runs after each step;
`verbose = true` prints the time, steps, last Δt, N and evaluations at each saved time (`progress_line(integ)` gives
that line when you drive the integrator yourself).
With [DotPlots](https://github.com/scottweady/DotPlots.jl) loaded, `dotplot(sol)` draws an `Evolution` (one frame
per saved time; `animate = true` plays them) and `dotplot(X)` a curve or a vector of curves (a package extension,
ext/ThetaLDotPlotsExt.jl).

## The problem, the integrator, the velocity protocol

```julia
function velocity(Xs, ev; warm)                      # all curves at once
    U = [κbar(Xs) .- curvature(X) for X in Xs]
    (; U, ctx = something_to_warm_start_the_next_call)   # or just return U
end
prob  = CurveProblem(Xs, velocity; symbol = q -> -q^2, rigid = :full, tags = [:a, :b])
integ = init(prob, ETDRK2())                          # evaluates the velocity at t0 (ev.kind == :init)
status = advance!(integ, 0.1; errtol = 1e-6)          # one step to t = 0.1 exactly: :success or :rejected
step!(integ, 0.01)                                    # one step of 0.01
sol = solve(prob, ETDRK2(); tspan = (0, 1), dt = 0.01, saveat = 0.1)   # or the driver
sol = solve(prob; tspan = (0, 1), adaptive = true, errtol = 1e-5)      # the step from the error estimate
```

- **`velocity(Xs, ev::Evaluation; warm)`** returns one vector per curve (the normal velocity at its grid
  points), or a named tuple `(; U, ctx)`; `ctx` comes back as `warm` at the next evaluation (a warm start). `ev`
  says what the evaluation is for: `t` (exact), `kind` (`:init`, `:stage`, `:end`; do expensive extras only at
  `:end`), `stage`, `c`, `h`, `attempt`, the curves' `ids` and `tags`, `topology` (bumped by every edit), and `base`
  (the states at the start of the step). Throw `StepFailure(reason)` to reject the step.
- **`symbol`** is the growth rate σ of a small-scale mode, θ_t ≈ σφ̂, as a function of the wavenumber along the
  curve q = 2π|k|/L: `q -> -q^2` (U = −κ), `q -> -S / 2 * q^3` (Hele-Shaw with surface tension S, scaled as in
  HLS 1994), `q -> -q^4` (U = κ_ss); or `(q, X) -> σ` when it depends on the curve X (anisotropy). σ < 0 damps;
  it is frozen per curve at the start of each step, and a positive σ gives a warning.
- **`rigid`**: the rigid motion (`:full`, `:translation`) taken out of U and integrated exactly; `:none` for the
  plain Hou–Lowengrub–Shelley step.
- **The integrator** (`ThetaLIntegrator`): `integ.t`, `integ.states`, `integ.U`, `integ.ctx`, `integ.ids`,
  `integ.tags`, `integ.rigid` (each curve's `RigidMotion` over the last step), `integ.stats` (steps, evaluations,
  rejections, `error_estimate`, `last_rejection`). Read them; change them only through the functions below.
  `integ.scheme`, `prob.velocity` and `prob.symbol` may be replaced between steps.
- **Edits**: `double_grid!(integ, i)`, `halve_grid!(integ, i)` (N changed, U kept), `set_curve!`, `add_curve!`
  (returns the new id), `remove_curve!`, `offset_curve!`, `curve_index(integ, id)`. Since F couples the curves,
  every other edit invalidates all velocities; the next `advance!` evaluates them first.
- **Rejection**: `advance!` returns `:rejected` and leaves the integrator exactly as it was when the velocity
  throws `StepFailure`, a stage has a non-positive length, the embedded error estimate exceeds `errtol` (checked
  before the end evaluation), or `prob.validate(states)` is false. A retry gets `ev.attempt = 2` and the committed
  context.
- **Restart**: `cp = checkpoint(integ)` is plain data (versioned; arrays, the ids and tags, the resolved filter);
  `init(prob, scheme; from = cp)` continues the run identically.
- **Stateless level**: `thetaL_step(scheme, prob, states, U, t, t_end; warm, …) -> StepResult`.
- **State carried with the curves**: `CurveProblem(…; bulk = BulkState(y0; symbol), fields = (; c =
  CurveField(c0; transport, symbol)))`. The velocity sees `ev.bulk` and `ev.fields.c` at each stage and returns
  `(; U, W, rates = (; c = …), ydot)`; the step advances y and c with the scheme (their stiff parts exactly, by
  their own symbols), rolls them back on rejection, carries the fields through the edits and checkpoints both.
  A field's `transport` along the curve (W the material tangential velocity, T the frame's): `:density`
  (c_t ∋ −((W − T)c)_s − c L_t/L, ∫c ds kept), `:material` (−(W − T)c_s) or `:none`. `integ.bulk`,
  `integ.fields.c`; `sol.bulk`, `sol.fields` from `solve`.

## Schemes

| Scheme | Order | Evaluations per step | Error estimate | With no symbol |
|:--|:--|:--|:--|:--|
| `ETDRK2()` (default) | 2 | 2 (midpoint, end) | embedded exponential Euler, before the end evaluation | explicit midpoint |
| `ExpEuler()` | 1 | 1 (end) | none | explicit Euler |

Both take `filter_order = 36` and `filter_cutoff` (`0`: N/2; an integer; `:auto`: from `noise_cutoff` after the
first evaluation). The end evaluation is the next step's first stage. Per curve and mode the step is

    stage i:  φ̂⁽ⁱ⁾ = e^{cᵢz} φ̂ + h Σⱼ Aᵢⱼ(z) Rⱼ,   end:  φ̂ₙ₊₁ = e^{z} φ̂ + h Σᵢ bᵢ(z) Rᵢ,   z = σh,

with R = θ_t − σφ̂ explicit, L by the same weights at z = 0, and x₀ by the averaged rigid motion (exact when it is
constant). Rigid motion is exact for every scheme (the rotation of a near-circle, which U barely shows, is left to
the deformation: `rigid_part`); the exponential weights keep a state at rest fixed under a symbol.

## Limits

- **Close approach.** The stiff part of F is assumed local to each curve and diagonal in k. Where two curves, or two
  parts of one, come within a gap d, the coupling is stiff at k ~ L/d in a way no symbol expresses, and the step
  needs Δt ≲ d·(scale). `min_distance(states)` monitors it.
- **Smoothness of F.** θ_t contains U_α by spectral differentiation: U must be smooth to `tol` at the grid points.
  A velocity interpolated from a finite element trace or a grid flux jump is not; fit it first.
- **Tangential transport.** After the rigid part is taken out, the tangential velocity of the deformation moves the
  points along the curve explicitly: Δt ≲ s_α/(k_max |T|). A first-order symbol damps it; steps beyond it are
  unstable (test/test_stepping.jl, the manufactured ellipse).
- **Frozen symbols** must bound the true stiffness (anisotropy: the maximum of γ + γ″ on the curve).
- **Initial layer.** From a shape far from equilibrium under a stiff symbol (surface diffusion, σ = −q⁴), the first steps
  are first order in Δt: the small scales decay within a step and drive L_t, which the explicit update of L follows
  only for Δt ≪ 1/|σ|. Second order once relaxed (test/test_stepping.jl); `adaptive = true` in `evolve`/`solve`
  lets the error estimate pick small steps there (`examples/step_control.jl`).

## Utilities

`min_distance(states)`, `curve_integral(X, f)`, `offset_curve(X, δ)` (the parallel curve), `remove_area_rate(states,
U; per_curve)` (per curve, or one constant for all: only the total area kept), `frame_velocity(X, U)` (the velocity
of the tracked points, for a mesh that follows them), `rigid_part(X, U)`, `noise_cutoff`, `spectral_tail`,
`double_grid`/`halve_grid`; geometry of a state: `curve_points`, `curve_function`, `tangent_angles`, `curvature`,
`arclength_derivative`, `area`, `centroid`, `perimeter`, `curve_moments`, `velocity_components`, `velocity_vectors`.

## Examples

`examples/` has its own environment: ThetaL by path, and [DotPlots](https://github.com/scottweady/DotPlots.jl)
from its repository, to draw the curves in the terminal (played in place, then overlaid):

| Example | Shows |
|:--|:--|
| `mean_curvature_flow.jl` | `evolve`, a star under U = −κ, `q -> -q^2`, Δt ∝ 1/κ² |
| `surface_diffusion.jl` | a flower relaxing under U = κ_ss (`arclength_derivative`), `q -> -q^4` |
| `rigid_rotation.jl` | a star turning about a point outside it, ten steps a turn, keeping its shape |
| `vortex_windup.jl` | a disk wound into a spiral by a vortex; N doubles as the curve sharpens |
| `coarsening.jl` | a coupled mean field, `init`/`advance!` by hand, particles removed as they vanish |
| `step_control.jl` | step size from the embedded estimate: `evolve(…; adaptive = true, errtol)`, and by hand (`advance!(…; errtol)`, `stats.error_estimate`) |
| `anisotropy.jl` | a curve-dependent symbol `(q, X) -> σ`, `solve` on a `CurveProblem` |
| `hele_shaw.jl` | an expanding bubble growing fingers: U from a boundary integral solve, `q -> -S * q^3`, N doubling |
| `vesicles.jl` | Willmore flow of membranes with area and length held: a pressure and a tension per curve, `q -> -q^4` |
| `hele_shaw_bubbles.jl` | three bubbles screening each other: one coupled velocity for all, Δt from `integ.U`, stopped at close approach |

```bash
julia --project=examples examples/coarsening.jl
```

## Extending

- **A scheme**: a struct `<: OneStepScheme` with one field `filter::SpectralFilter`, a keyword constructor, and
  methods `stage_nodes(scheme)` (the cᵢ, c₁ = 0), `scheme_weights(scheme, z) -> (A, b, b̂)` (tuples of the weights at
  z = σh) and optionally `has_error_estimate(scheme) = true` (src/stepping/schemes.jl). The generic step does the rest.

## Layout

| Path | |
|:--|:--|
| `src/ThetaL.jl` | module docstring (equations, conventions, code map), exports, includes, the PrecompileTools workload |
| `src/grid/` | `fourier.jl` (periodic grid functions: wavenumbers, spectral derivative and antiderivative, `double_grid`, `halve_grid`, a trigonometric interpolant), `work.jl` (their work arrays: in-place FFT plans and scratch vectors per N, kept per task) |
| `src/curves/` | `state.jl` (the state and its geometry), `construct.jl` (`from_curve`, `from_arclength_samples`), `measures.jl` (utilities for callers) |
| `src/stepping/` | `symbols.jl` (the rates σ(q) of a curve), `filter.jl` (`SpectralFilter`, `noise_cutoff`), `kinematics.jl` (`thetaL_rhs`, `rigid_part`, `frame_velocity`), `schemes.jl` (`ETDRK2`, `ExpEuler`, their weights), `carried.jl` (the step's part for carried state: transport, combination), `step.jl` (the generic step, `thetaL_step`) |
| `src/problem/` | `evaluation.jl` (`Evaluation`, `StepFailure`, reading the velocity's result), `problem.jl` (`CurveProblem`), `carried.jl` (`CurveField`, `BulkState`: state carried with the curves) |
| `src/integrator/` | `integrator.jl` (`ThetaLIntegrator`, `init`, `advance!`, `step!`), `edits.jl`, `checkpoint.jl` |
| `src/drivers/` | `evolve.jl` (`evolve`, `solve`, `Evolution`, `save_curves`) |
| `ext/ThetaLDotPlotsExt.jl` | `dotplot` of a curve, of curves and of an `Evolution`, loaded with DotPlots (a weak dependency) |
| `examples/` | the examples above, and their environment |
| `test/` | `test_curves.jl`, `test_stepping.jl` (kinematics and schemes), `test_integrator.jl`, `test_drivers.jl`, `test_carried.jl` (state carried with the curves); `common/` (shared helpers) |
| `docs/` | the documentation: `make.jl`, `src/` (pages, assets), `figures/` (`make_data.jl`, `plot_figures.py`, the data), `codemap.jl` (the code map) |

Compile time: the kernels are loops on `Vector{Float64}`/`Vector{ComplexF64}` (no fused broadcasts or
comprehensions, each of which compiles on first use), and the functions a caller passes in are left unspecialized
(`@nospecialize`), so one compiled step per scheme serves every curve, velocity and symbol.

## Tests

From this directory (one file: `test_args=["integrator"]`):

```bash
julia -t 8 --project -e 'using Pkg; Pkg.test()'
```

Dependencies: FFTW, StaticArrays, PrecompileTools, and the standard library LinearAlgebra. Not registered: use it by
path, `Pkg.develop(path = "path/to/ThetaL")`.
