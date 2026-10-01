# CIS Project Part 1 - Inverted Pendulum

This project designs a discrete-time LQR controller and an angle-only state
observer for the nonlinear actuated inverted pendulum in the
[Part 1 assignment](CIS_Project_Part1-1.pdf). The main entry point is
[`pendulum_core.m`](pendulum_core.m). It contains the model, controller,
observer, required experiments, metrics, assertions, and plots in one file.

## Running the project

1. Open MATLAB and make this repository the current folder.
2. Run:

   ```matlab
   pendulum_core
   ```

3. Inspect the generated `results/` folder and the summary printed in the
   Command Window.

The four result figures also remain open as normal MATLAB windows. Use the
figure tabs to switch between them, or run `close all` when you are finished.

The equation-based workflow runs with base MATLAB. If Control System Toolbox
is present, the script uses `dlqr` and `place`; otherwise it uses the included
two-state fallback calculations.

The four main figures are:

- `results/01_open_loop.png`
- `results/02_controller_experiments.png`
- `results/03_observer_experiments.png`
- `results/04_final_disturbed_angle_only.png`

The two main evidence tables are:

- `results/controller_metrics.csv`
- `results/observer_metrics.csv`

`results/all_metrics.csv` is a convenience export containing all experiment
rows in one table.

### Optional Simulink/Simscape model

After installing Simulink, Simscape, and Simscape Multibody, run:

```matlab
build_pendulum_model
pendulum_core
```

The builder creates `actuated_pendulum_final.slx` without changing the supplied
team models. When that model is available, the last section of
`pendulum_core.m` attempts the same final angle-only disturbed case in
Simulink. A successful validation additionally creates:

- `results/05_simulink_validation.png`
- `results/simulink_metrics.csv`

The equation-based MATLAB workflow has been run. The optional final
Simulink/Simscape validation is a separate check and is **not yet claimed as
tested** here.

## Model and controller equations

The nonlinear physical model is

\[
J\ddot{\theta}=mgl\sin(\theta)-b\dot{\theta}+u,
\qquad J=ml^2,
\]

where `theta = 0` is upright. With
`x = [theta; omega]` and `omega = theta_dot`, MATLAB integrates

\[
\dot{x}_1=x_2,
\qquad
\dot{x}_2=\frac{mgl\sin(x_1)-bx_2+u}{J}.
\]

Near upright, `sin(theta) approximately theta`, giving

\[
\dot{x}=A_cx+B_cu,
\quad
A_c=\begin{bmatrix}0&1\\mgl/J&-b/J\end{bmatrix},
\quad
B_c=\begin{bmatrix}0\\1/J\end{bmatrix}.
\]

The exact zero-order-hold discretization is computed from

\[
\exp\left(
\begin{bmatrix}A_c&B_c\\0&0\end{bmatrix}T_s
\right)
=
\begin{bmatrix}A_d&B_d\\0&1\end{bmatrix}.
\]

The full-state LQR uses `uCommand = -K*x`. The angle-only version uses
`uCommand = -K*xhat`. Actuator saturation and the unknown plant disturbance
are applied as

\[
u_{act}=\operatorname{sat}(u_{command},u_{min},u_{max}),
\qquad
u_{plant}=u_{act}+v.
\]

Only angle is measured:

\[
y[k]=Cx[k]+w[k], \qquad C=\begin{bmatrix}1&0\end{bmatrix}.
\]

The predictor observer is

\[
\hat{x}[k+1]=A_d\hat{x}[k]+B_du_{act}[k]
+L\left(y[k]-C\hat{x}[k]\right).
\]

### Why ZOH and Unit Delay are needed

The physical pendulum is continuous, but the digital controller updates only
at `t_k = k*Ts`. A Zero-Order Hold keeps the selected torque constant between
updates, exactly matching the assumption used to calculate `Ad` and `Bd`.

The observer needs its previous estimate `xhat[k]` to calculate
`xhat[k+1]`. In Simulink, the Unit Delay stores that estimate for one sample.
It is the observer's discrete memory and also prevents an algebraic loop.

## Nominal design and numerical anchors

The nominal values are:

| Quantity | Value |
|---|---:|
| Mass `m` | `1 kg` |
| Pivot-to-COM distance `l` | `0.5 m` |
| Damping `b` | `0.05 N*m*s/rad` |
| Moment of inertia `J` | `0.25 kg*m^2` |
| Sampling period `Ts` | `0.02 s` |
| Torque limits | `-10` to `+10 N*m` |
| LQR weights | `Q = diag([100 1])`, `R = 1` |
| Observer poles | `0.60`, `0.55` |

Expected regression anchors are:

- Continuous open-loop poles: approximately `+4.3306` and `-4.5306`.
- Discrete LQR gain: approximately `K = [14.8022  2.76935]`.
- Controller poles: approximately `0.88531 +/- 0.05158i`.
- Observer gain: approximately `L = [0.853845; 9.39513]`.
- The nominal full-state `35 deg` case settles in about `0.74 s`, reaches
  about `9.04 N*m`, and does not saturate.
- The nominal angle-only case starts from the deliberately wrong estimate
  `xhat0 = [0;0]`, settles in about `0.8 s`, and briefly reaches the
  `10 N*m` limit while the estimate converges.

Small numerical differences are normal across MATLAB releases and solver
settings. The assertions use physical tolerances rather than exact equality.

## What each script section does

1. **Physical plant and nominal design model** - defines separate `p` (the
   actual nonlinear plant) and `m` (the controller design model), plus `Ts`,
   weights, limits, observer poles, and initial estimate.
2. **Nonlinear open loop** - runs `ode45` at the exact upright state and from
   `+1 deg` and `-1 deg` perturbations. It saves the time histories and the
   eye-shaped phase portrait in `01_open_loop.png`.
3. **Linearization, ZOH, LQR, and observer** - forms `Ac`, `Bc`, `Ad`, `Bd`,
   `K`, and `L`; verifies controllability, observability, and pole stability;
   and prints the nominal matrices and poles.
4. **Common experiment definition** - creates one reproducible base
   configuration used by all later cases.
5. **Q/R trade-off** - compares nominal weights, a larger angle weight, and a
   larger control penalty using full-state feedback from `10 deg`.
6. **Initial conditions and actuator limits** - compares `5`, `35`, and
   `120 deg` with a `10 N*m` limit, plus a deliberately unsuccessful
   `35 deg` case limited to `2 N*m`.
7. **Observer poles** - compares slow `[0.85 0.80]`, nominal `[0.60 0.55]`,
   and fast `[0.25 0.20]` observers, all starting with an incorrect estimate.
8. **Sampling time** - redesigns both controller and observer at `0.01`,
   `0.02`, and `0.10 s`; it compares full-state and observer feedback while
   keeping approximately the same continuous observer speed.
9. **Measurement noise and finite difference** - adds seeded `0.5 deg`
   standard-deviation angle noise and compares the model-based observer with
   `(y[k]-y[k-1])/Ts` as a velocity estimate.
10. **Model uncertainty** - designs the controller using `0.8` times the real
    mass, `1.2` times the real length, and `1.5` times the real damping, then
    applies it to the unchanged physical plant.
11. **Input disturbance/final demonstration** - runs angle-only feedback with
    a seeded `0.20 N*m` standard-deviation unknown torque disturbance added
    after saturation.
12. **Metrics tables** - writes the two main CSV summaries and the combined
    convenience CSV.
13. **Controller figure** - saves Q/R, torque, initial-condition/limit, and
    model-mismatch comparisons in `02_controller_experiments.png`.
14. **Observer figure** - saves observer-error, sampling-time, and noisy
    finite-difference comparisons in `03_observer_experiments.png`.
15. **Final angle-only figure** - saves true/estimated angle and velocity plus
    actuator, disturbance, and total plant torque in
    `04_final_disturbed_angle_only.png`.
16. **Optional Simulink verification** - runs `actuated_pendulum_final.slx`
    when it exists and the required products are available; otherwise it
    skips this check without preventing the MATLAB experiments.

The local functions after Section 16 implement parameter construction,
experiment execution, nonlinear sampled-data simulation, metric calculation,
LQR/observer fallbacks, logged-signal extraction, and optional Simulink plots.

## How to interpret the required experiments

### Open-loop stability

The exact state `[0;0]` remains still because its derivative is exactly zero.
That does **not** mean the equilibrium is stable. The positive pole and the
growth of either small perturbation show that upright is unstable. Damping
eventually produces loops around the hanging configurations near `+/-180 deg`.

### Q and R

Increasing the angle weight produces faster correction but requests more
torque and uses more control effort. Increasing `R` makes torque more
expensive, giving a gentler but slower response. These trends are not absolute
once saturation is active: two aggressive designs can become nearly identical
while both are clipped.

### Initial state and actuator limits

The nominal `10 N*m` controller stabilizes all selected initial angles,
including the tested `120 deg` case. The `2 N*m` case does not settle and
instead reaches the hanging region. This demonstrates an empirical region of
attraction; it does not prove that `10 N*m` can recover every possible state.
The static comparison with `m*g*l` is useful intuition but is not a global
recoverability theorem.

### Input disturbance

The disturbance is unknown to the observer and is added after actuator
saturation. The final case remains close to upright, but angle and torque have
small persistent fluctuations. Report last-second RMS as well as settling
time, because a continuously disturbed response is not expected to become
identically zero.

### Model mismatch

Only the controller/observer design model is changed; the nonlinear physical
plant remains nominal. The tested moderate combined mismatch still settles,
with a slightly changed transient and effort. That supports robustness to this
specific mismatch only; it is not proof of robustness to arbitrary parameter
errors.

### Observer poles and incorrect initial estimate

All three clean, noise-free observer choices converge. Faster discrete poles
(smaller magnitudes) remove initial estimation error sooner, while slower poles
take longer. Very fast observers generally transmit more measurement noise
into the state estimate and torque, so the fastest possible poles are not
automatically best. The clean pole sweep and the separate noisy experiment
should be interpreted together rather than claiming one universally optimal
choice.

### Sampling time

The controller and observer are redesigned at each tested `Ts`, so this is a
fair sampled-data comparison rather than reusing gains designed at `0.02 s`.
All three selected sampling periods stabilize in the current run, but the
`0.10 s` observer case has a larger excursion and slower response. In general,
a longer sampling period lets the unstable plant move farther between updates;
the three tested points do not establish a universal failure threshold.

### Measurement noise

Angle noise produces nonzero estimation and control fluctuations even after
the initial transient. Use the matched post-transient errors and last-second
RMS values when comparing methods; a whole-run RMSE can be dominated by the
deliberately incorrect initial estimate.

### Finite-difference velocity

Finite difference is simple and reacts without a plant model, but it divides
the difference of two noisy samples by `Ts`, so smaller sampling periods can
strongly amplify noise. The observer uses the model and applied actuator
command to filter information over time, but its result depends on model
accuracy and pole selection. Therefore the comparison is a trade-off, not a
claim that either method is always superior.

## Meaning of the main metrics

- **Settling time:** first time after which both `|theta| <= 1 deg` and
  `|omega| <= 0.1 rad/s` remain satisfied.
- **Settled:** whether that permanent entry occurs within the simulation.
- **Maximum angle:** largest absolute raw angle reached.
- **Maximum torque:** largest absolute saturated actuator command.
- **Control effort:** `Ts*sum(uActuator.^2)`.
- **Saturation percentage:** fraction of commands that requested torque beyond
  the actuator limits.
- **Estimation RMSE:** angle/velocity disagreement between `x` and `xhat`.
- **Last-second RMS:** steady behavior used especially for noise and
  disturbance cases.

The script asserts mathematical invariants and a few intended reference
outcomes. A deliberately failing experimental case is recorded as data rather
than treated as a software error.
