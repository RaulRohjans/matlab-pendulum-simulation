# CIS Project Part 2 - Cart-pole swing-up and stabilization

This folder contains a complete implementation of the Part 2 brief:

1. Exact zero-order-hold discretization of the supplied upright linear model.
2. Saturated discrete LQR stabilization of the nonlinear plant.
3. An empirical four-state capture-zone study.
4. Force-limited iLQR swing-up with backward recursion, nonlinear forward
   rollout, backtracking line search, convergence checks, and final gains.
5. Sampled tracking with a latched switch to upright LQR.
6. Nonlinear Simscape validation, feasibility checks, figures, and metrics.

The original supplied files are unchanged in `template/`. Siddhi's independent
attempt is preserved in `attempts/Siddhi/`. Keep those directories off the
MATLAB path because their filenames intentionally overlap the canonical files.

## Run the complete workflow

Make this `part2` directory the MATLAB Current Folder, then run:

```matlab
result = setup_cartpole;
```

This starts at the downward equilibrium (`theta = 180 deg`), determines the
capture criterion, optimizes the swing-up, checks the nominal trajectory, runs
`cartpole.slx`, switches to LQR after capture, and opens the result figures.

`cartpole.slx` is generated from the untouched supplied model with MATLAB's
supported Stateflow API. If it ever needs to be regenerated, run:

```matlab
build_cartpole_model
```

## Run the Simulink model directly

Open `cartpole.slx` and click **Run**. The model contains a complete validated
default swing-up design in its model workspace, so no setup command is needed
first. Multibody Explorer opens automatically and shows the moving cart and
pendulum in 3-D. The `x` and `u` signals are also returned in the simulation
output.

Running `setup_cartpole` remains the recommended complete workflow because it
recomputes the controller, verifies feasibility, opens the plots, and
temporarily overrides the embedded defaults without changing the model.

For numerical-only or headless execution, disable the 3-D window explicitly:

```matlab
result = setup_cartpole(struct('showMultibodyExplorer',false));
```

For an equation-only check without Simulink:

```matlab
result = setup_cartpole(struct( ...
    'runSimulink',false, ...
    'makePlots',false));
```

If Java-based Multibody Explorer is unavailable, numerical Simscape validation
still works. The driver disables that optional viewer. A simple MATLAB
animation is also available after a successful run:

```matlab
animate_cartpole(result)
```

## Required experiments

Reproduce the nonlinear capture-zone study:

```matlab
capture = capture_zone_experiment;
```

Compare force and cart limits, weights, horizons, and zero versus sinusoidal
initial guesses:

```matlab
experiments = run_part2_experiments;
```

The experiment script reports optimization failures instead of stopping at the
first unsuccessful case. Note that `xmax` is a feasibility check in the
provided plant, not a physical hard stop; the cart-position cost shapes the
trajectory to respect it.

## Main files

- `setup_cartpole.m` - complete LQR, iLQR, feasibility, Simscape, and plotting
  workflow.
- `backward_lqr.m` - supplied LQ tracking recursion.
- `forward_pass.m` - nonlinear rollout with force clipping before propagation.
- `trajectory_cost.m` - objective from the project brief.
- `ilqr.m` - line search, stopping criteria, regularization, and final gains.
- `linearize_trajectory.m` - central-difference discrete Jacobians.
- `controller_block.m` - readable copy of the controller embedded in the SLX.
- `capture_zone_experiment.m` - empirical saturated-LQR capture criterion.
- `run_part2_experiments.m` - required sensitivity experiments.

## Angle convention

Optimization keeps the nominal angle unwrapped so the quadratic objective and
linearization remain smooth. The real tracking error, upright LQR state, and
capture check wrap angles using

```matlab
atan2(sin(angle),cos(angle))
```

so configurations separated by an integer multiple of `2*pi` are treated as
the same physical pole orientation. If the nominal trajectory ends without
entering the capture zone, the embedded controller stops tracking and applies
zero force; post-processing reports the failure explicitly.

## Generated evidence

Outputs are written to `results/` and remain visible as MATLAB figure windows:

- `01_capture_zone.png`
- `02_ilqr_nominal.png`
- `03_simscape_validation.png`
- `04_ilqr_experiments.png`
- `capture_zone_summary.csv`
- `part2_metrics.csv`
- `ilqr_experiments.csv`
