# API reference

```@docs
ThetaL
```

## Curves

```@docs
ThetaLState
from_curve
from_arclength_samples
curve_function
curve_points
tangent_angles
curvature
arclength_derivative
velocity_components
velocity_vectors
curve_moments
area
centroid
perimeter
spectral_tail
double_grid
halve_grid
```

## Utilities

```@docs
min_distance
curve_integral
offset_curve
remove_area_rate
frame_velocity
rigid_part
RigidMotion
rigid_velocity
```

## The problem and the velocity protocol

```@docs
CurveProblem
CurveField
BulkState
Evaluation
StepFailure
```

## Schemes and the step

```@docs
OneStepScheme
ETDRK2
ExpEuler
SpectralFilter
noise_cutoff
thetaL_rhs
thetaL_step
StepResult
```

For scheme authors (not exported; see [Adding a scheme](@ref)):

```@docs
ThetaL.stage_nodes
ThetaL.scheme_weights
ThetaL.has_error_estimate
ThetaL.with_filter
ThetaL.evaluations_per_step
```

## The integrator

```@docs
ThetaLIntegrator
ThetaL.IntegratorStats
init(::CurveProblem, ::OneStepScheme)
advance!
step!(::ThetaLIntegrator, ::Real)
progress_line
double_grid!
halve_grid!
set_curve!
add_curve!
remove_curve!
offset_curve!
curve_index
checkpoint
```

## Drivers

```@docs
evolve
solve(::CurveProblem, ::OneStepScheme)
Evolution
save_curves
```
