# Control and Intelligent Systems project

The repository is split by assignment part so that the two physical models,
controllers, templates, and result sets cannot shadow one another on the
MATLAB path.

- [`part1/`](part1/) - actuated inverted-pendulum modelling, LQR, observer,
  experiments, and the final Simulink/Simscape model.
- [`part2/`](part2/) - underactuated cart-pole swing-up with iLQR followed by
  saturated upright LQR stabilization.

The root contains only repository-wide files. Do not use `addpath(genpath(...))`
because both parts deliberately retain untouched templates and team attempts
with duplicate function names. Instead, make the desired part the MATLAB
Current Folder before running it.

## Run Part 1

```matlab
cd('<repository>/part1')
pendulum_core
```

## Run Part 2

```matlab
cd('<repository>/part2')
setup_cartpole
```

See each part's README for its assumptions, experiments, and outputs.
