# Getting started

## Installation

ThetaL is not registered (Julia 1.10 or later; dependencies FFTW, StaticArrays,
PrecompileTools, and the standard library LinearAlgebra). Use it from another environment by path:

```julia
using Pkg
Pkg.develop(path = "path/to/ThetaL")
```

The examples have an environment of their own, which also brings [DotPlots](https://github.com/scottweady/DotPlots.jl)
to draw the curves in the terminal:

```
cd ThetaL
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples examples/mean_curvature_flow.jl
```

## One call: `evolve`

```julia
using ThetaL

X = ξ -> (1 + 0.3cos(5ξ)) .* (cos(ξ), sin(ξ))     # any counterclockwise 2π-periodic parametrization
sol = evolve(X, (X, t) -> -curvature(X), (0, 0.1);
             dt = 1e-3, symbol = q -> -q^2, saveat = 0.02)

sol.t                      # [0, 0.02, …, 0.1]: the saved times, landed on exactly
sol(0.1)                   # the curve then, a ThetaLState
area(sol(0.1)), perimeter(sol(0.1)), centroid(sol(0.1))
save_curves("star.csv", sol)                          # t, curve id, x, y
```

With [DotPlots](https://github.com/scottweady/DotPlots.jl) loaded (`using DotPlots`), `dotplot(sol)` draws the saved
curves in the terminal, overlaid and colored by time; `dotplot(sol; animate = true)` plays them first, and
`dotplot(sol(0.1))` draws one time. The keywords are DotPlots' (`width`, `every`, …).

- The curve may be a parametrization, a `ThetaLState`, or a vector of them (several curves, with `U(Xs, t)`
  returning one velocity per curve).
- `U(X, t)` returns the outward normal velocity at the grid points of the curve `X`, computed from its geometry:
  `curve_points(X)`, `curvature(X)`, `tangent_angles(X)`, `arclength_derivative(X, f)`. For a velocity field
  ``u(x)``, `velocity_components(X, u.(curve_points(X)))[1]` is ``u \cdot n``.
- `symbol` is the stiff part of `U` (below); `q -> -q^2` is that of ``U = -\kappa``.
- `dt` is a number or a function `dt(X, t)`. N is chosen to the tolerance `tol` and doubled when the curve needs
  it; `coarsen = true` also halves it once the curve has become smooth.
- `verbose = true` prints a line at each saved time: the time, the steps, the last Δt, N and the velocity
  evaluations, with the wall time; and a line whenever N changes. `progress_line(integ)` is that line, for a loop of your own.

## A problem and its integrator

For anything more (an expensive velocity, curves that come and go, your own control of the step), build a
`CurveProblem` and drive the integrator yourself:

```julia
function velocity(Xs, ev; warm)                           # the velocities of all curves at once
    κbar = 2π * length(Xs) / sum(perimeter, Xs)
    [κbar .- curvature(X) for X in Xs]                    # or (; U, ctx) to hand a context to the next call
end
prob  = CurveProblem(Xs, velocity; symbol = q -> -q^2)
integ = init(prob, ETDRK2())                              # evaluates the velocity at t = 0

while integ.t < 1
    status = advance!(integ, integ.t + 0.01)              # one step, landing exactly; :success or :rejected
    status === :success || break
end
integ.states, integ.U, integ.ids, integ.stats
```

`solve(prob, ETDRK2(); tspan, dt, saveat)` runs the same loop as `evolve` on a problem.

## The velocity protocol

`velocity(Xs, ev; warm)` is called with all curves, and told what the evaluation is for:

| `ev.` | |
|:--|:--|
| `t` | the time; at the end of a step, exactly the time stepped to |
| `kind` | `:init` (the start, or after a curve edit), `:stage` (inside a step), `:end` (the end of a step) |
| `stage`, `c`, `h` | the stage, its fraction of the step, the step |
| `attempt` | 1, or more when the step is retried after a rejection |
| `ids`, `tags` | per curve, a stable id and your own tag (ThetaL never reads it) |
| `topology` | bumped by every curve edit and change of N: a warm start keyed on it knows when to start over |
| `base` | the states at the start of the step |

It returns one vector per curve, or a named tuple `(; U, ctx)`: `ctx` comes back as `warm` at the next evaluation.
An expensive solve warm-starts from it, and does extra work (an interior field, diagnostics) only at
`ev.kind === :end`. To reject a step (a solver that did not converge), throw `StepFailure(reason)`.

## The stiff symbol

The symbol is the growth rate ``\sigma`` of a small-scale mode of the curve, ``\theta_t \approx \sigma\hat\phi``,
as a function of its wavenumber along the curve, ``q = 2\pi|k|/L``: `symbol = q -> -q^2` for curve shortening,
`q -> -q^4` for surface diffusion, `q -> -S / 2 * q^3` for Hele-Shaw with surface tension ``S`` (in the scaling of
Hou, Lowengrub & Shelley 1994). ``\sigma < 0``
damps the mode. When it depends on the curve, give `symbol = (q, X) -> …`, with `X` the curve (a `ThetaLState`); it
is frozen per curve at the start of each step. A table or a quadrature works as well as a formula (see
[Symbols, filter and resolution](@ref)).

## State carried with the curves

A problem may carry state that moves with the curves: a **bulk state** (one vector y off the curves, such as the
unknowns of a field with memory) and **fields on the curves** (a value per grid point of each curve, such as a
concentration c). The step advances them with the curves and with the same scheme, the stiff part of each treated
exactly by its own symbol. A rejected step leaves them as they were; checkpoints keep them.

```julia
prob = CurveProblem(curves, velocity; symbol = q -> -q^2,
                    bulk = BulkState(y0; symbol = σy),                         # σy: growth rates, diagonal
                    fields = (; c = CurveField(c0; transport = :density, symbol = q -> -D * q^2)))

function velocity(Xs, ev; warm)
    y = ev.bulk                                  # the bulk state at this stage
    c = ev.fields.c                              # the field at this stage, one vector per curve
    …                                            # the solve
    (; U, W, rates = (; c = S), ydot = ẏ)        # W: the material tangential velocity, one vector per curve
end
integ = init(prob); advance!(integ, 0.1)
integ.bulk, integ.fields.c                       # at integ.t; `solve` keeps them at the saved times: sol.bulk, sol.fields
```

A field slides along the grid where the material moves tangentially (at W) differently from the grid points (at
the frame's T, which keeps them at equal arclength). ThetaL adds that transport:

| `transport` | the step adds to the field's rate S |
|:--|:--|
| `:density` (default) | a conserved density per unit length: ``-((W - T)c)_s - c\,L_t/L``; ``\int c\,ds`` is kept exactly on the grid |
| `:material` | a value carried by the material points: ``-(W - T)\,c_s`` |
| `:none` | nothing: S is the whole ``c_t`` at the grid points |

so `rates` holds only the physics (a source, adsorption, ``D c_{ss}`` with its symbol). Under a rigid motion
W = T and a field is exactly unchanged. The fields follow their curve through the edits: interpolated by
`double_grid!` and `halve_grid!`, resampled by `set_curve!`, given to `add_curve!` (`fields = (; c = …)`, or the
field's `init` when that is a number or a function `X -> values`), dropped by `remove_curve!`. `evolve` stays the
one-call driver for curves alone; carried state goes through `CurveProblem` and `solve` or the integrator.

## Editing curves, rejection, restarts

```julia
id = add_curve!(integ, ξ -> (3 + 0.5cos(ξ), 0.5sin(ξ)); tag = :new)   # returns the new id
remove_curve!(integ, curve_index(integ, id))
double_grid!(integ, 1); halve_grid!(integ, 1)                       # change N; the velocities are kept
offset_curve!(integ, 1, 0.01); set_curve!(integ, 1, new_state)

advance!(integ, t; errtol = 1e-6)     # rejected if the embedded error estimate exceeds errtol
integ.stats.last_rejection            # :error_estimate, :nonpositive_length, :invalid_state, or your reason

cp = checkpoint(integ)                # plain data (versioned)
integ = init(prob, ETDRK2(); from = cp)
```

Because the velocity of every curve may depend on every curve, any edit other than a change of N makes the next
`advance!` evaluate all velocities first. A rejected step leaves the integrator exactly as it was; the retry gets
`ev.attempt = 2` and the committed context.
