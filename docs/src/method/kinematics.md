# Kinematics and the rigid split

## A curve in θ–L variables

A closed counterclockwise curve ``\Gamma`` of length ``L`` is parametrized by ``\alpha \in [0, 2\pi)`` at equal
arclength, ``s_\alpha = L/2\pi``. Its tangent angle is ``\theta(\alpha) = \alpha + \phi(\alpha)`` with ``\phi``
periodic, and the curve is

```math
x(\alpha) = x_0 + \frac{L}{2\pi} \int_0^\alpha (\cos\theta, \sin\theta) \, d\alpha'.
```

The state of a curve (`ThetaLState`) is ``\phi`` at ``N`` equispaced ``\alpha``, the length ``L`` and the point
``x_0 = x(0)``. The tangent is ``\tau = (\cos\theta, \sin\theta)``, the outward normal
``n = (\sin\theta, -\cos\theta)``, and the curvature ``\kappa = \theta_s = (1 + \phi_\alpha)/s_\alpha``, positive on
convex parts. A state is made from any parametrization by `from_curve` (spectrally, with N chosen to a
tolerance) or from points already at equal arclength by `from_arclength_samples`.

## The motion

The points move with ``x_t = U n + T \tau``: ``U`` is the given normal velocity, and the tangential velocity
``T`` is chosen to keep the points at equal arclength (Hou, Lowengrub & Shelley 1994):

```math
L_t = \int_0^{2\pi} U \, \theta_\alpha \, d\alpha,
```

```math
T_\alpha = \frac{L_t}{2\pi} - U \theta_\alpha, \qquad T(0) = 0,
```

```math
\theta_t = \frac{T \theta_\alpha - U_\alpha}{s_\alpha}, \qquad x_{0,t} = U(0)\, n(0).
```

`thetaL_rhs(X, U, filt)` computes ``\theta_t``, ``L_t`` and ``x_{0,t}``.

## The rigid split

With ``T(0) = 0`` the points slide along a curve that turns, at a rate set by the tangential velocity at ``x_0``,
which changes as ``x_0`` slides; and ``x_0`` itself moves on a circle that an explicit step cannot follow exactly.
So no time integrator of these equations keeps a rigid rotation exact. ThetaL splits the normal velocity of each
curve, at every stage of every step, into its rigid part and the rest:

```math
U = U_{\mathrm{rig}} + U_d, \qquad U_{\mathrm{rig}} = V \cdot n - \Omega\, (x - x_c) \cdot \tau,
```

where the translation ``V`` and rotation rate ``\Omega`` about the center ``x_c`` (the mean of the grid points) give
the rigid normal velocity closest to ``U`` on the grid (least squares; `rigid_part`). The rigid part moves the
curve without moving the points along it: it adds ``\Omega`` to ``\theta_t`` and carries ``x_0`` on the rigid
motion, ``z' = i\Omega z + a`` in complex form, which the step integrates exactly. The deformation ``U_d`` goes
through the equations above; it moves ``x_0`` by ``d = U_d(0)\, n(0)``, a velocity of the curve that turns with
it, so each stage's ``d`` is taken back to the frame of the start of the step (``e^{-i\Theta} d``, ``\Theta`` the
rotation that stage's state has made), combined with the scheme's weights, and turned with the step's rotation.
Under a rigid motion ``\hat\phi_k`` (``k \neq 0``) and ``L`` stay fixed and only the mean of ``\phi`` moves, by
``\Omega h``, which every scheme integrates exactly: translation and rotation are exact to roundoff at any time
step ([Method overview](@ref), panel (a)).

The split needs only ``U``, not the material tangential velocity, since the tangential motion does not change the
curve. A curve whose shape is steady while material flows along it (an edge current) has ``U = 0`` and stays
exactly steady. Near a circle the rotation hardly shows in ``U`` (a circle turning about its center has ``U = 0``),
so noise in ``U`` is fitted as a rotation rate of any size, different at every stage. That is harmless: the rotation
of a near-circle about its center hardly moves it, and as ``x_0``'s drift turns with the curve, the points only
slide along it (test/test_stepping.jl, a near-circle with noise in ``U``). When the problem transports fields along
the curves, the velocity gives the material tangential velocity ``W``, and the rigid part is fitted to the whole
velocity ``u = U n + W\tau`` instead: ``V`` the mean of ``u`` and
``\Omega = \sum (x - x_c) \times u / \sum |x - x_c|^2``, well determined for every curve, so the points follow the
material. `CurveProblem(...; rigid = :translation)` takes out the translation only; `rigid = :none` gives the plain
Hou–Lowengrub–Shelley step.

## The velocity of the tracked points

A bulk mesh that follows the grid points needs their velocity, not only its normal part: `frame_velocity(X, U)`
returns ``U n + T\tau`` with the tangential velocity of the step (the rigid motion's, plus the equal-arclength one of
the deformation).
