"""
    ThetaL

Boundary evolution: θ–L motion of closed plane curves (Hou, Lowengrub & Shelley 1994). N curves move with outward
normal velocities that depend on all of them, (U₁, …, U_N) = F(X₁, …, X_N; t), given by the caller as one function;
the stiff part of each curve's motion, whose rate per Fourier mode (the small-scale symbol) the caller gives too, is
treated exactly by an exponential integrator, and the rigid part of each curve's motion is integrated exactly.
Nothing here knows how the velocity is computed.

A curve Γ (counterclockwise) is represented by its tangent angle θ(α) and length L, with α ∈ [0, 2π) an
equal-arclength parameter (s_α = L/2π), sampled at N equispaced α. The state holds the periodic part φ = θ − α,
the length L and the point x₀ = x(α = 0). A state is made from any parametrization of the curve (`from_curve`,
spectral) or from points already at equal arclength (`from_arclength_samples`).

Conventions: τ = (cos θ, sin θ), outward normal n = (sin θ, −cos θ), κ = θ_s. The points move with
x_t = U n + T τ (U the given normal velocity, T chosen to keep equal arclength):
  L_t = ∫₀^{2π} U θ_α dα,    T_α = L_t/2π − U θ_α,  T(0) = 0,
  θ_t = (T θ_α − U_α)/s_α,   x₀,t = U(0) n(0).
Rigid motion: the normal velocity is split as U = U_rig + U_d, with U_rig = V·n − Ω (x − x_c)·τ the rigid motion
closest to U (`rigid_part`). The tangential velocity is the rigid one plus that of U_d from the equation above, so a
rigid motion does not move the points along the curve: it adds Ω to θ_t, and moves x₀ on the rigid motion, which the
step integrates exactly. A rigid translation and rotation are therefore exact at any Δt.
Small scales: θ_t ≈ σ(q) φ̂_k with σ the symbol, the growth rate of a mode of wavenumber q = 2π|k|/L (σ < 0
damps), frozen per curve at the start of each step. The one-step schemes (`ETDRK2`, `ExpEuler`) are exponential
Runge–Kutta methods: the remainder R = θ_t − σφ̂ is explicit, integrated against e^{σt} with ϕ-function weights.
A spectral filter is applied to φ after each stage and to U_d before differentiation.

Code map (src/):
- `grid/`: periodic grid functions (`fourier.jl`: wavenumbers, spectral derivative, `double_grid`, `halve_grid`;
  real FFTs on half spectra, complex ones for curves as x + iy) and their work arrays (`work.jl`: FFT plans and
  scratch vectors per N, kept per task).
- Each state caches what is computed from it more than once (its half spectrum, θ_α, its points); the step makes
  its states with their spectra, so a step needs five real FFTs per stage, and one complex pair per new state.
- `curves/`: the state (`state.jl`), states from curves (`construct.jl`), utilities for callers (`measures.jl`).
- `stepping/`: the stiff symbol (`symbols.jl`), the spectral filter (`filter.jl`), the kinematics and the rigid
  split (`kinematics.jl`), the schemes and their weights (`schemes.jl`), the generic step (`step.jl`).
- `problem/`: the velocity protocol (`evaluation.jl`: `Evaluation`, `StepFailure`), `CurveProblem`, and the state
  carried with the curves (`carried.jl`: `CurveField`, `BulkState`; its step in `stepping/carried.jl`).
- `integrator/`: `ThetaLIntegrator`, `init`, `advance!`, `step!` (`integrator.jl`), curve edits (`edits.jl`),
  `checkpoint` (`checkpoint.jl`).
- `drivers/`: `evolve`, `solve`, `Evolution`, `save_curves`.
Outside src/: `ext/ThetaLDotPlotsExt.jl`, `dotplot` of curves and of an `Evolution` when DotPlots is loaded.
"""
module ThetaL

using StaticArrays
using LinearAlgebra: dot, mul!
using PrecompileTools: @setup_workload, @compile_workload
import FFTW

# curves: states and their geometry
export ThetaLState, from_curve, from_arclength_samples, curve_function, curve_points, tangent_angles, curvature,
       arclength_derivative, velocity_components, velocity_vectors, curve_moments, area, centroid, perimeter,
       spectral_tail, double_grid, halve_grid
# utilities for callers
export min_distance, curve_integral, offset_curve, remove_area_rate, frame_velocity
# the problem and the velocity protocol
export CurveProblem, CurveField, BulkState, Evaluation, StepFailure, RigidMotion, rigid_velocity, rigid_part
# schemes and the step
export OneStepScheme, ETDRK2, ExpEuler, SpectralFilter, noise_cutoff, thetaL_rhs, thetaL_step, StepResult
# the integrator and its edits
export ThetaLIntegrator, init, advance!, step!, progress_line, double_grid!, halve_grid!, set_curve!, add_curve!,
       remove_curve!, offset_curve!, curve_index, checkpoint
# drivers
export evolve, solve, Evolution, save_curves

const Vec2 = SVector{2,Float64}

# periodic grid functions and their work arrays
include("grid/fourier.jl")
include("grid/work.jl")
# the state of a curve, and states from a curve
include("curves/state.jl")
include("curves/construct.jl")
include("curves/measures.jl")
# the step: symbol, filter, kinematics, schemes, the generic step
include("stepping/symbols.jl")
include("stepping/filter.jl")
include("stepping/kinematics.jl")
include("stepping/schemes.jl")
include("problem/carried.jl")
include("problem/evaluation.jl")
include("problem/problem.jl")
include("stepping/step.jl")
include("stepping/carried.jl")
# the integrator
include("integrator/integrator.jl")
include("integrator/edits.jl")
include("integrator/checkpoint.jl")
# the drivers
include("drivers/evolve.jl")

# What nearly every script runs, run while the package precompiles (PrecompileTools) so that the package image
# holds what it compiles, including what is reached by runtime dispatch. The caller's functions are left
# unspecialized, so this one compiled chain serves every curve, velocity and symbol of a script. Kept to that chain:
# each further call here costs precompile time after every change to ThetaL, while a function left out compiles
# at its first use, in the scripts that use it (the curve edits, checkpoints, ExpEuler, the utilities).
@setup_workload begin
    circle = ξ -> (cos(ξ), sin(ξ))
    @compile_workload begin
        X = from_curve(circle, 64)
        curvature(X); curve_points(X); area(X); velocity_components(X, curve_points(X))
        sol = evolve(circle, (X, t) -> -curvature(X), (0, 0.02); dt=0.01, symbol=q -> -q^2, saveat=0.01)
        sol(0.01)
        integ = init(CurveProblem([circle], (Xs, ev; warm) -> [-curvature(X) for X in Xs]; symbol=q -> -q^2))
        advance!(integ, 0.01); progress_line(integ)
    end
end

end
