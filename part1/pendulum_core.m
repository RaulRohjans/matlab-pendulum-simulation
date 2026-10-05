%% Actuated inverted pendulum - complete Project Part 1 workflow
% This single script covers the nonlinear model, linearization, exact ZOH
% discretization, full-state LQR, angle-only observer, and the experiments
% requested in the Part 1 brief.  All random experiments use fixed seeds.
%
% Outputs are written to this Part 1 folder's results subfolder:
%   all_metrics.csv, controller_metrics.csv, observer_metrics.csv
%   01_open_loop.png ... 04_final_disturbed_angle_only.png
%
% The equation-based work runs with base MATLAB. Control System Toolbox is
% used when available; small two-state fallbacks are included for dlqr and
% place. If actuated_pendulum_final.slx exists here, it is run at the end.

clear; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(projectRoot, 'results');
if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

rng(7, 'twister');
g = 9.81;

%% 1. Physical plant and nominal design model
% p is the actual nonlinear plant. m is the model used for controller and
% observer design. Keeping them separate makes the uncertainty test clear.
p = makeParameters(1.0, 0.5, 0.05);
p.theta0 = 35;                  % degrees, for the Simulink model
p.omega0 = 0;                   % rad/s
m = makeParameters(p.m, p.l, p.b);

Ts = 0.02;
Tend = 6;
umin = -10;
umax = 10;
Q = diag([100, 1]);
R = 1;
C = [1, 0];
xhat0 = [0; 0];
observerPoles = [0.60, 0.55];
runSimulinkModel = true;
% Figures are also exported to results/, but remain open in small MATLAB
% windows so the response can be inspected immediately after the run.

odeOptions = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%% 2. Nonlinear open-loop model: equilibrium and small perturbations
openLoopTend = 5;
freeDynamics = @(time, state) pendulumDynamics( ...
    time, state, 0, p, g);

[tEquilibrium, xEquilibrium] = ode45( ...
    freeDynamics, [0, openLoopTend], [0; 0], odeOptions);
[tPositive, xPositive] = ode45( ...
    freeDynamics, [0, openLoopTend], [deg2rad(1); 0], odeOptions);
[tNegative, xNegative] = ode45( ...
    freeDynamics, [0, openLoopTend], [-deg2rad(1); 0], odeOptions);

assert(max(abs(xEquilibrium(:))) < 1e-12, ...
    'The exact upright equilibrium should remain at rest when u = 0.');
assert(max(abs(xPositive(:,1))) > deg2rad(30), ...
    'The positive perturbation did not demonstrate open-loop instability.');
assert(max(abs(xNegative(:,1))) > deg2rad(30), ...
    'The negative perturbation did not demonstrate open-loop instability.');

fig = figure(1); clf(fig);
set(fig, 'Name', '1 - Open-loop pendulum', 'NumberTitle', 'off', ...
    'Visible', 'on', 'Color', 'w', 'Position', [100, 100, 900, 650]);
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(layout);
plot(ax, tEquilibrium, rad2deg(xEquilibrium(:,1)), ...
    'Color', [0.25, 0.25, 0.25], 'LineWidth', 1.4, ...
    'DisplayName', 'Exact equilibrium');
hold(ax, 'on');
plot(ax, tPositive, rad2deg(xPositive(:,1)), ...
    'LineWidth', 1.4, 'DisplayName', '+1 deg perturbation');
plot(ax, tNegative, rad2deg(xNegative(:,1)), ...
    'LineWidth', 1.4, 'DisplayName', '-1 deg perturbation');
yline(ax, 180, 'k:', 'HandleVisibility', 'off');
yline(ax, -180, 'k:', 'HandleVisibility', 'off');
grid(ax, 'on');
xlabel(ax, 'Time [s]');
ylabel(ax, 'Angle from upright [deg]');
title(ax, 'Nonlinear pendulum without control');
legend(ax, 'Location', 'best');

ax = nexttile(layout);
plot(ax, rad2deg(xPositive(:,1)), xPositive(:,2), ...
    'LineWidth', 1.4, 'DisplayName', '+1 deg perturbation');
hold(ax, 'on');
plot(ax, rad2deg(xNegative(:,1)), xNegative(:,2), ...
    'LineWidth', 1.4, 'DisplayName', '-1 deg perturbation');
grid(ax, 'on');
xlabel(ax, 'Angle from upright [deg]');
ylabel(ax, 'Angular velocity [rad/s]');
title(ax, 'Phase portrait');
legend(ax, 'Location', 'best');

exportgraphics(fig, fullfile(resultsDir, '01_open_loop.png'), ...
    'Resolution', 160);

%% 3. Linearization, stability, exact ZOH model, LQR and observer
% With sin(theta) approximately theta near the upright equilibrium:
%   xdot = Ac*x + Bc*u,
%   Ac = [0 1; m*g*l/J -b/J], Bc = [0; 1/J].
[Ac, Bc] = continuousLinearModel(m, g);
continuousPoles = eig(Ac);
assert(any(real(continuousPoles) > 0), ...
    'The upright linearized model should have one unstable pole.');

nominalDesign = designDiscreteSystem( ...
    m, g, Ts, Q, R, C, observerPoles);
Ad = nominalDesign.Ad;
Bd = nominalDesign.Bd;
K = nominalDesign.K;
L = nominalDesign.L;

fprintf('\nContinuous-time open-loop poles:\n');
disp(continuousPoles);
fprintf('Nominal Ad:\n'); disp(Ad);
fprintf('Nominal Bd:\n'); disp(Bd);
fprintf('Nominal LQR gain K:\n'); disp(K);
fprintf('Nominal controller poles:\n'); disp(nominalDesign.controllerPoles);
fprintf('Nominal observer gain L:\n'); disp(L);
fprintf('Nominal observer poles:\n'); disp(nominalDesign.observerPoles);

%% 4. Common experiment definition
baseConfig = makeBaseConfig(m, Ts, Tend, Q, R, umin, umax, ...
    observerPoles, xhat0);

%% 5. Q/R trade-off experiment (full state, small initial angle)
qConfigs = repmat(baseConfig, 1, 3);
for index = 1:3
    qConfigs(index).mode = 'full-state';
    qConfigs(index).group = 'Q and R';
    qConfigs(index).x0 = [deg2rad(10); 0];
end
qConfigs(1).name = 'Q/R nominal';
qConfigs(2).name = 'larger angle weight';
qConfigs(2).Q = diag([400, 1]);
qConfigs(3).name = 'larger control penalty';
qConfigs(3).R = 10;
qResults = runCases(p, g, qConfigs, odeOptions);

assert(all(arrayfun(@(item) item.metrics.Settled, qResults)), ...
    'At least one Q/R comparison case did not stabilize.');

%% 6. Initial-condition and actuator-limit experiment
initialConfigs = repmat(baseConfig, 1, 4);
initialAngles = [5, 35, 120, 35];
limits = [10, 10, 10, 2];
initialNames = {'5 deg, limit 10', '35 deg, limit 10', ...
    '120 deg, limit 10', '35 deg, limit 2'};
for index = 1:4
    initialConfigs(index).name = initialNames{index};
    initialConfigs(index).group = 'initial state and limits';
    initialConfigs(index).mode = 'full-state';
    initialConfigs(index).x0 = [deg2rad(initialAngles(index)); 0];
    initialConfigs(index).umin = -limits(index);
    initialConfigs(index).umax = limits(index);
end
initialResults = runCases(p, g, initialConfigs, odeOptions);

assert(initialResults(2).metrics.Settled, ...
    'The nominal full-state case did not stabilize from 35 degrees.');
assert(~initialResults(4).metrics.Settled, ...
    ['The deliberately tight 2 N*m limit unexpectedly met the settling ' ...
     'criterion; choose a smaller limit to preserve the contrast.']);

%% 7. Observer-pole experiment with an incorrect initial estimate
poleConfigs = repmat(baseConfig, 1, 3);
poleSets = {[0.85, 0.80], [0.60, 0.55], [0.25, 0.20]};
poleNames = {'slow observer', 'nominal observer', 'fast observer'};
for index = 1:3
    poleConfigs(index).name = poleNames{index};
    poleConfigs(index).group = 'observer poles';
    poleConfigs(index).mode = 'observer';
    poleConfigs(index).x0 = [deg2rad(35); 0];
    poleConfigs(index).xhat0 = [0; 0];
    poleConfigs(index).observerPoles = poleSets{index};
end
poleResults = runCases(p, g, poleConfigs, odeOptions);
cleanObserverResult = poleResults(2);

assert(cleanObserverResult.metrics.Settled, ...
    'The nominal angle-only observer case did not stabilize.');
assert(isfinite(cleanObserverResult.metrics.EstimationSettlingTime_s), ...
    'The nominal state estimate did not converge.');

%% 8. Sampling-time experiment: full-state and observer feedback
sampleTimes = [0.01, 0.02, 0.10];
samplingConfigs = repmat(baseConfig, 1, 2*numel(sampleTimes));
writeIndex = 0;
for timeIndex = 1:numel(sampleTimes)
    localTs = sampleTimes(timeIndex);
    % Keep approximately the same continuous observer speed at each Ts.
    mappedObserverPoles = exp([-30, -25]*localTs);
    for modeIndex = 1:2
        writeIndex = writeIndex + 1;
        samplingConfigs(writeIndex).Ts = localTs;
        samplingConfigs(writeIndex).observerPoles = mappedObserverPoles;
        samplingConfigs(writeIndex).x0 = [deg2rad(35); 0];
        samplingConfigs(writeIndex).group = 'sampling time';
        if modeIndex == 1
            samplingConfigs(writeIndex).mode = 'full-state';
            samplingConfigs(writeIndex).name = sprintf( ...
                'Ts %.3f full state', localTs);
        else
            samplingConfigs(writeIndex).mode = 'observer';
            samplingConfigs(writeIndex).name = sprintf( ...
                'Ts %.3f observer', localTs);
        end
    end
end
samplingResults = runCases(p, g, samplingConfigs, odeOptions);

%% 9. Measurement-noise experiment and alternative velocity estimate
measurementNoiseConfig = baseConfig;
measurementNoiseConfig.name = 'angle noise 0.5 deg';
measurementNoiseConfig.group = 'measurement noise';
measurementNoiseConfig.mode = 'observer';
measurementNoiseConfig.x0 = [deg2rad(35); 0];
measurementNoiseConfig.measurementNoiseStd = deg2rad(0.5);
measurementNoiseConfig.seed = 202;
measurementNoiseResult = runCases( ...
    p, g, measurementNoiseConfig, odeOptions);

assert(measurementNoiseResult.metrics.AngleRMSLast1s_deg < 3, ...
    'The noisy-measurement case has excessive steady-state angle error.');
assert(measurementNoiseResult.metrics.OmegaEstimateRMSE_rad_s ...
        < measurementNoiseResult.metrics.FiniteDifferenceOmegaRMSE_rad_s, ...
    ['After the two-second transient, the observer should estimate velocity ' ...
     'more accurately than an unfiltered finite difference.']);

fprintf(['Finite-difference velocity RMSE under angle noise: %.4f rad/s\n' ...
         'Observer velocity RMSE under angle noise:          %.4f rad/s\n'], ...
    measurementNoiseResult.metrics.FiniteDifferenceOmegaRMSE_rad_s, ...
    measurementNoiseResult.metrics.OmegaEstimateRMSE_rad_s);

%% 10. Model-uncertainty experiment
mismatchModel = makeParameters(0.8*p.m, 1.2*p.l, 1.5*p.b);
uncertaintyConfig = baseConfig;
uncertaintyConfig.name = 'mismatched design model';
uncertaintyConfig.group = 'model uncertainty';
uncertaintyConfig.mode = 'full-state';
uncertaintyConfig.model = mismatchModel;
uncertaintyConfig.x0 = [deg2rad(35); 0];
uncertaintyResult = runCases(p, g, uncertaintyConfig, odeOptions);

assert(uncertaintyResult.metrics.Settled, ...
    'The selected moderate model mismatch did not stabilize.');

%% 11. Input-disturbance experiment and final angle-only demonstration
disturbanceConfig = baseConfig;
disturbanceConfig.name = 'angle only with input disturbance';
disturbanceConfig.group = 'input disturbance';
disturbanceConfig.mode = 'observer';
disturbanceConfig.x0 = [deg2rad(35); 0];
disturbanceConfig.inputDisturbanceStd = 0.20;
disturbanceConfig.seed = 101;
disturbanceResult = runCases(p, g, disturbanceConfig, odeOptions);

assert(disturbanceResult.metrics.AngleRMSLast1s_deg < 2, ...
    'The disturbed angle-only case has excessive steady-state angle error.');

%% 12. Save compact metrics tables
allResults = [qResults, initialResults, poleResults, samplingResults, ...
    measurementNoiseResult, uncertaintyResult, disturbanceResult];
controllerResults = [qResults, initialResults, uncertaintyResult, ...
    disturbanceResult];
observerResults = [poleResults, samplingResults, ...
    measurementNoiseResult, disturbanceResult];

allMetrics = resultsToTable(allResults);
controllerMetrics = resultsToTable(controllerResults);
observerMetrics = resultsToTable(observerResults);

writetable(allMetrics, fullfile(resultsDir, 'all_metrics.csv'));
writetable(controllerMetrics, ...
    fullfile(resultsDir, 'controller_metrics.csv'));
writetable(observerMetrics, ...
    fullfile(resultsDir, 'observer_metrics.csv'));

fprintf('\nController experiment summary:\n');
disp(controllerMetrics(:, {'Case', 'SettlingTime_s', 'Settled', ...
    'MaxAngle_deg', 'MaxTorque_Nm', 'ControlEffort_Nm2s', ...
    'MaxPlantTorque_Nm', 'PlantInputEffort_Nm2s'}));
fprintf('\nObserver experiment summary:\n');
disp(observerMetrics(:, {'Case', 'Ts_s', 'Settled', ...
    'AngleEstimateRMSE_deg', 'OmegaEstimateRMSE_rad_s'}));

%% 13. Controller experiment figure
colors = lines(4);
fig = figure(2); clf(fig);
set(fig, 'Name', '2 - Controller experiments', 'NumberTitle', 'off', ...
    'Visible', 'on', 'Color', 'w', 'Position', [100, 100, 1100, 750]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(layout);
hold(ax, 'on');
for index = 1:numel(qResults)
    plot(ax, qResults(index).t, rad2deg(qResults(index).x(1,:)), ...
        'LineWidth', 1.3, 'Color', colors(index,:), ...
        'DisplayName', qResults(index).config.name);
end
grid(ax, 'on');
xlabel(ax, 'Time [s]'); ylabel(ax, 'Angle [deg]');
title(ax, 'Effect of Q and R'); legend(ax, 'Location', 'best');

ax = nexttile(layout);
hold(ax, 'on');
for index = 1:numel(qResults)
    stairs(ax, qResults(index).t(1:end-1), ...
        qResults(index).uActuator, 'LineWidth', 1.2, ...
        'Color', colors(index,:), ...
        'DisplayName', qResults(index).config.name);
end
grid(ax, 'on');
xlabel(ax, 'Time [s]'); ylabel(ax, 'Actuator torque [N*m]');
title(ax, 'Q/R control effort'); legend(ax, 'Location', 'best');

ax = nexttile(layout);
hold(ax, 'on');
for index = 1:numel(initialResults)
    plot(ax, initialResults(index).t, ...
        rad2deg(initialResults(index).x(1,:)), ...
        'LineWidth', 1.2, 'Color', colors(index,:), ...
        'DisplayName', initialResults(index).config.name);
end
grid(ax, 'on');
xlabel(ax, 'Time [s]'); ylabel(ax, 'Angle [deg]');
title(ax, 'Initial states and actuator limits');
legend(ax, 'Location', 'best');

ax = nexttile(layout);
plot(ax, initialResults(2).t, ...
    rad2deg(initialResults(2).x(1,:)), 'LineWidth', 1.3, ...
    'DisplayName', 'Nominal model');
hold(ax, 'on');
plot(ax, uncertaintyResult.t, rad2deg(uncertaintyResult.x(1,:)), ...
    '--', 'LineWidth', 1.3, 'DisplayName', 'Mismatched design model');
grid(ax, 'on');
xlabel(ax, 'Time [s]'); ylabel(ax, 'Angle [deg]');
title(ax, 'Model uncertainty'); legend(ax, 'Location', 'best');

exportgraphics(fig, fullfile(resultsDir, ...
    '02_controller_experiments.png'), 'Resolution', 160);

%% 14. Observer experiment figure
fig = figure(3); clf(fig);
set(fig, 'Name', '3 - Observer experiments', 'NumberTitle', 'off', ...
    'Visible', 'on', 'Color', 'w', 'Position', [100, 100, 1100, 750]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(layout);
hold(ax, 'on');
for index = 1:numel(poleResults)
    plot(ax, poleResults(index).t, rad2deg( ...
        poleResults(index).x(1,:) - poleResults(index).xhat(1,:)), ...
        'LineWidth', 1.2, 'Color', colors(index,:), ...
        'DisplayName', poleResults(index).config.name);
end
grid(ax, 'on');
xlabel(ax, 'Time [s]'); ylabel(ax, 'Angle error [deg]');
title(ax, 'Observer-pole comparison'); legend(ax, 'Location', 'best');

ax = nexttile(layout);
hold(ax, 'on');
for index = 1:numel(poleResults)
    plot(ax, poleResults(index).t, ...
        poleResults(index).x(2,:) - poleResults(index).xhat(2,:), ...
        'LineWidth', 1.2, 'Color', colors(index,:), ...
        'DisplayName', poleResults(index).config.name);
end
grid(ax, 'on');
xlabel(ax, 'Time [s]'); ylabel(ax, 'Velocity error [rad/s]');
title(ax, 'Velocity-estimation convergence');
legend(ax, 'Location', 'best');

ax = nexttile(layout);
hold(ax, 'on');
observerSampleIndices = find(arrayfun(@(item) ...
    strcmp(item.config.mode, 'observer'), samplingResults));
for index = 1:numel(observerSampleIndices)
    resultIndex = observerSampleIndices(index);
    plot(ax, samplingResults(resultIndex).t, ...
        rad2deg(samplingResults(resultIndex).x(1,:)), ...
        'LineWidth', 1.2, 'Color', colors(index,:), ...
        'DisplayName', samplingResults(resultIndex).config.name);
end
grid(ax, 'on');
xlabel(ax, 'Time [s]'); ylabel(ax, 'Angle [deg]');
title(ax, 'Sampling-time comparison'); legend(ax, 'Location', 'best');

ax = nexttile(layout);
plot(ax, measurementNoiseResult.t, ...
    measurementNoiseResult.x(2,:) - measurementNoiseResult.xhat(2,:), ...
    'LineWidth', 1.2, 'DisplayName', 'Observer error');
hold(ax, 'on');
plot(ax, measurementNoiseResult.t, ...
    measurementNoiseResult.x(2,:) ...
    - measurementNoiseResult.omegaFiniteDifference, ...
    'LineWidth', 1.0, 'DisplayName', 'Finite-difference error');
grid(ax, 'on');
xlim(ax, [2, measurementNoiseResult.config.Tend]);
xlabel(ax, 'Time [s]'); ylabel(ax, 'Velocity error [rad/s]');
title(ax, 'Noisy angle: observer vs finite difference');
legend(ax, 'Location', 'best');

exportgraphics(fig, fullfile(resultsDir, ...
    '03_observer_experiments.png'), 'Resolution', 160);

%% 15. Final disturbed angle-only result figure
fig = figure(4); clf(fig);
set(fig, 'Name', '4 - Final disturbed angle-only response', ...
    'NumberTitle', 'off', 'Visible', 'on', 'Color', 'w', ...
    'Position', [100, 100, 900, 750]);
layout = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(layout);
plot(ax, disturbanceResult.t, rad2deg(disturbanceResult.x(1,:)), ...
    'LineWidth', 1.4, 'DisplayName', 'True');
hold(ax, 'on');
plot(ax, disturbanceResult.t, rad2deg(disturbanceResult.xhat(1,:)), ...
    '--', 'LineWidth', 1.2, 'DisplayName', 'Estimated');
grid(ax, 'on'); ylabel(ax, 'Angle [deg]');
title(ax, 'Angle-only feedback with plant-input disturbance');
legend(ax, 'Location', 'best');

ax = nexttile(layout);
plot(ax, disturbanceResult.t, disturbanceResult.x(2,:), ...
    'LineWidth', 1.4, 'DisplayName', 'True');
hold(ax, 'on');
plot(ax, disturbanceResult.t, disturbanceResult.xhat(2,:), ...
    '--', 'LineWidth', 1.2, 'DisplayName', 'Estimated');
grid(ax, 'on'); ylabel(ax, 'Velocity [rad/s]');
legend(ax, 'Location', 'best');

ax = nexttile(layout);
stairs(ax, disturbanceResult.t(1:end-1), ...
    disturbanceResult.uActuator, 'LineWidth', 1.2, ...
    'DisplayName', 'Actuator torque');
hold(ax, 'on');
stairs(ax, disturbanceResult.t(1:end-1), ...
    disturbanceResult.inputDisturbance, ':', 'LineWidth', 1.0, ...
    'DisplayName', 'External disturbance');
stairs(ax, disturbanceResult.t(1:end-1), ...
    disturbanceResult.uPlant, '--', 'LineWidth', 1.0, ...
    'DisplayName', 'Total plant torque');
yline(ax, disturbanceConfig.umax, 'k:', 'HandleVisibility', 'off');
yline(ax, disturbanceConfig.umin, 'k:', 'HandleVisibility', 'off');
grid(ax, 'on'); xlabel(ax, 'Time [s]'); ylabel(ax, 'Torque [N*m]');
legend(ax, 'Location', 'best');

exportgraphics(fig, fullfile(resultsDir, ...
    '04_final_disturbed_angle_only.png'), 'Resolution', 160);

%% 16. Optional final Simulink/Simscape verification
% A Zero-Order Hold samples the measurements/controller output at t_k and
% keeps u[k] constant on [t_k,t_(k+1)). The observer Unit Delay stores
% xhat[k] while the MATLAB Function block computes xhat[k+1]; it provides
% discrete memory and prevents an algebraic loop.
%
% The model is optional so all mathematical experiments remain runnable on
% machines without Simulink or Simscape Multibody.
finalModelFile = fullfile(projectRoot, 'actuated_pendulum_final.slx');
if runSimulinkModel && isfile(finalModelFile)
    % Newer MATLAB releases can report Simulink's `sim` entry point as a
    % package/folder (exist == 7), so only treat it as unavailable when it
    % cannot be resolved at all.
    if isempty(which('sim'))
        warning('Simulink is unavailable; the final SLX check was skipped.');
    else
        try
            [~, finalModelName] = fileparts(finalModelFile);
            load_system(finalModelFile);

            % Use the same final disturbed angle-only case as above.
            reference = disturbanceResult;
            m = reference.config.model;
            Ts = reference.config.Ts;
            Tend = reference.config.Tend;
            Q = reference.config.Q;
            R = reference.config.R;
            umin = reference.config.umin;
            umax = reference.config.umax;
            xhat0 = reference.config.xhat0;
            observerPoles = reference.config.observerPoles;
            Ad = reference.design.Ad;
            Bd = reference.design.Bd;
            C = reference.design.C;
            K = reference.design.K;
            L = reference.design.L;
            % The final model uses uniform random sources with half-width
            % amplitudes. A uniform signal on [-a,a] has std = a/sqrt(3).
            disturbanceAmplitude = sqrt(3) ...
                * reference.config.inputDisturbanceStd;
            measurementNoiseAmplitude = sqrt(3) ...
                * reference.config.measurementNoiseStd;
            useObserver = true;

            simulationInput = Simulink.SimulationInput(finalModelName);
            % These names exactly match build_pendulum_model.m.
            variableNames = {'p', 'Ts', 'Tend', 'K', 'umin', 'umax', ...
                'Ad', 'Bd', 'C', 'L', 'xhat0', 'useObserver', ...
                'disturbanceAmplitude', 'measurementNoiseAmplitude'};
            variableValues = {p, Ts, Tend, K, umin, umax, ...
                Ad, Bd, C, L, xhat0, useObserver, ...
                disturbanceAmplitude, measurementNoiseAmplitude};
            for variableIndex = 1:numel(variableNames)
                simulationInput = simulationInput.setVariable( ...
                    variableNames{variableIndex}, ...
                    variableValues{variableIndex});
            end
            simulationInput = simulationInput.setModelParameter( ...
                'StopTime', num2str(Tend), ...
                'ReturnWorkspaceOutputs', 'on');
            simulationOutput = sim(simulationInput);
            fprintf('Final Simulink model ran successfully: %s\n', ...
                finalModelFile);

            [thetaTime, thetaData] = extractLoggedSignal( ...
                simulationOutput, {'theta_sim', 'theta', 'theta_true'});
            [omegaTime, omegaData] = extractLoggedSignal( ...
                simulationOutput, {'omega_sim', 'omega', 'omega_true'});
            [thetaHatTime, thetaHatData] = extractLoggedSignal( ...
                simulationOutput, {'theta_hat', 'theta_estimate'});
            [omegaHatTime, omegaHatData] = extractLoggedSignal( ...
                simulationOutput, {'omega_hat', 'omega_estimate'});
            [actuatorTorqueTime, actuatorTorqueData] = extractLoggedSignal( ...
                simulationOutput, {'u_sim', 'u_applied', 'torque', 'u'});
            [plantTorqueTime, plantTorqueData] = extractLoggedSignal( ...
                simulationOutput, {'u_plant_sim', 'u_plant', ...
                'plant_torque'});

            if ~isempty(thetaData)
                assert(abs(thetaData(end)) < deg2rad(5), ...
                    'The final Simulink angle did not finish near upright.');

                finalOmega = NaN;
                if ~isempty(omegaData)
                    finalOmega = omegaData(end);
                end
                maxActuatorTorque = NaN;
                actuatorEffort = NaN;
                if ~isempty(actuatorTorqueData)
                    maxActuatorTorque = max(abs(actuatorTorqueData));
                    actuatorEffort = loggedSignalEffort( ...
                        actuatorTorqueTime, actuatorTorqueData);
                end
                maxPlantTorque = NaN;
                plantInputEffort = NaN;
                if ~isempty(plantTorqueData)
                    maxPlantTorque = max(abs(plantTorqueData));
                    plantInputEffort = loggedSignalEffort( ...
                        plantTorqueTime, plantTorqueData);
                end
                simulinkMetrics = table(rad2deg(thetaData(end)), ...
                    finalOmega, max(abs(rad2deg(thetaData))), ...
                    maxActuatorTorque, actuatorEffort, maxPlantTorque, ...
                    plantInputEffort, ...
                    'VariableNames', {'FinalAngle_deg', ...
                    'FinalOmega_rad_s', 'MaxAngle_deg', ...
                    'MaxActuatorTorque_Nm', 'ActuatorEffort_Nm2s', ...
                    'MaxPlantTorque_Nm', 'PlantInputEffort_Nm2s'});
                writetable(simulinkMetrics, fullfile(resultsDir, ...
                    'simulink_metrics.csv'));

                saveSimulinkFigure(resultsDir, thetaTime, thetaData, ...
                    thetaHatTime, thetaHatData, omegaTime, omegaData, ...
                    omegaHatTime, omegaHatData, actuatorTorqueTime, ...
                    actuatorTorqueData, plantTorqueTime, plantTorqueData);
            else
                warning(['The final model ran but no recognized angle log ' ...
                    'was found. Expected theta_sim, theta, or theta_true.']);
            end
        catch simulinkError
            warning('Final Simulink verification failed: %s', ...
                simulinkError.message);
        end
    end
elseif runSimulinkModel
    fprintf(['Optional final model not found, so the Simulink check was ' ...
        'skipped:\n  %s\n'], finalModelFile);
end

fprintf('\nPart 1 equation-based workflow completed. Results are in:\n  %s\n', ...
    resultsDir);

%% Local functions
function parameters = makeParameters(mass, lengthToCOM, damping)
parameters.m = mass;
parameters.l = lengthToCOM;
parameters.b = damping;
parameters.J = mass*lengthToCOM^2;
end

function config = makeBaseConfig(model, Ts, Tend, Q, R, umin, umax, ...
        observerPoles, xhat0)
config.name = 'base case';
config.group = 'base';
config.mode = 'observer';
config.model = model;
config.Ts = Ts;
config.Tend = Tend;
config.Q = Q;
config.R = R;
config.umin = umin;
config.umax = umax;
config.x0 = [deg2rad(35); 0];
config.xhat0 = xhat0;
config.observerPoles = observerPoles;
config.inputDisturbanceStd = 0;
config.measurementNoiseStd = 0;
config.seed = 7;
end

function results = runCases(plant, gravity, configs, odeOptions)
numberOfCases = numel(configs);
results = repmat(emptyResult(), 1, numberOfCases);
for caseIndex = 1:numberOfCases
    results(caseIndex) = simulateClosedLoop( ...
        plant, gravity, configs(caseIndex), odeOptions);
end
end

function result = emptyResult()
result.config = struct();
result.plant = struct();
result.design = struct();
result.t = [];
result.x = [];
result.xhat = [];
result.y = [];
result.measurementNoise = [];
result.uCommand = [];
result.uActuator = [];
result.inputDisturbance = [];
result.uPlant = [];
result.omegaFiniteDifference = [];
result.metrics = struct();
end

function result = simulateClosedLoop(plant, gravity, config, odeOptions)
design = designDiscreteSystem(config.model, gravity, config.Ts, ...
    config.Q, config.R, [1, 0], config.observerPoles);

numberOfSteps = round(config.Tend/config.Ts);
t = (0:numberOfSteps)*config.Ts;
numberOfSamples = numel(t);

x = zeros(2, numberOfSamples);
xhat = zeros(2, numberOfSamples);
y = zeros(1, numberOfSamples);
uCommand = zeros(1, numberOfSteps);
uActuator = zeros(1, numberOfSteps);
uPlant = zeros(1, numberOfSteps);
x(:,1) = config.x0;
xhat(:,1) = config.xhat0;

previousRandomState = rng;
rng(config.seed, 'twister');
measurementNoise = config.measurementNoiseStd*randn(1, numberOfSamples);
inputDisturbance = config.inputDisturbanceStd*randn(1, numberOfSteps);
rng(previousRandomState);

for sampleIndex = 1:numberOfSteps
    y(sampleIndex) = design.C*x(:,sampleIndex) ...
        + measurementNoise(sampleIndex);

    if strcmp(config.mode, 'full-state')
        feedbackState = x(:,sampleIndex);
    elseif strcmp(config.mode, 'observer')
        feedbackState = xhat(:,sampleIndex);
    else
        error('Unknown feedback mode: %s', config.mode);
    end

    uCommand(sampleIndex) = -design.K*feedbackState;
    uActuator(sampleIndex) = min(max( ...
        uCommand(sampleIndex), config.umin), config.umax);

    % The external torque disturbance acts after actuator saturation. The
    % observer knows the actuator command but not the disturbance.
    uPlant(sampleIndex) = uActuator(sampleIndex) ...
        + inputDisturbance(sampleIndex);

    innovation = y(sampleIndex) - design.C*xhat(:,sampleIndex);
    xhat(:,sampleIndex + 1) = design.Ad*xhat(:,sampleIndex) ...
        + design.Bd*uActuator(sampleIndex) + design.L*innovation;

    heldInputDynamics = @(time, state) pendulumDynamics( ...
        time, state, uPlant(sampleIndex), plant, gravity);
    [~, stateDuringSample] = ode45(heldInputDynamics, ...
        [t(sampleIndex), t(sampleIndex + 1)], ...
        x(:,sampleIndex), odeOptions);
    x(:,sampleIndex + 1) = stateDuringSample(end,:)';
end

y(end) = design.C*x(:,end) + measurementNoise(end);
omegaFiniteDifference = [0, diff(y)/config.Ts];
metrics = calculateMetrics(config, plant, design, t, x, xhat, ...
    uCommand, uActuator, uPlant, omegaFiniteDifference);

result = emptyResult();
result.config = config;
result.plant = plant;
result.design = design;
result.t = t;
result.x = x;
result.xhat = xhat;
result.y = y;
result.measurementNoise = measurementNoise;
result.uCommand = uCommand;
result.uActuator = uActuator;
result.inputDisturbance = inputDisturbance;
result.uPlant = uPlant;
result.omegaFiniteDifference = omegaFiniteDifference;
result.metrics = metrics;
end

function metrics = calculateMetrics(config, plant, design, t, x, xhat, ...
        uCommand, uActuator, uPlant, omegaFiniteDifference)
angleTolerance = deg2rad(1);
velocityTolerance = 0.1;
insideSettlingBand = abs(x(1,:)) <= angleTolerance ...
    & abs(x(2,:)) <= velocityTolerance;
settlingTime = firstPermanentEntry(t, insideSettlingBand);

angleEstimationError = x(1,:) - xhat(1,:);
velocityEstimationError = x(2,:) - xhat(2,:);
insideEstimationBand = abs(angleEstimationError) <= deg2rad(0.5) ...
    & abs(velocityEstimationError) <= 0.1;
estimationSettlingTime = firstPermanentEntry(t, insideEstimationBand);

lastSecondState = t >= max(0, t(end) - 1);
lastSecondInput = t(1:end-1) >= max(0, t(end) - 1);
postTransientState = t >= 2;
saturated = uCommand < config.umin | uCommand > config.umax;

metrics.SettlingTime_s = settlingTime;
metrics.Settled = isfinite(settlingTime);
metrics.FinalAngle_deg = rad2deg(x(1,end));
metrics.FinalOmega_rad_s = x(2,end);
metrics.MaxAngle_deg = max(abs(rad2deg(x(1,:))));
metrics.MaxTorque_Nm = max(abs(uActuator));
metrics.ControlEffort_Nm2s = config.Ts*sum(uActuator.^2);
metrics.MaxPlantTorque_Nm = max(abs(uPlant));
metrics.PlantInputEffort_Nm2s = config.Ts*sum(uPlant.^2);
metrics.Saturation_pct = 100*mean(saturated);
metrics.AngleEstimateRMSE_deg = sqrt(mean( ...
    rad2deg(angleEstimationError).^2));
metrics.OmegaEstimateRMSE_rad_s = sqrt(mean( ...
    velocityEstimationError(postTransientState).^2));
metrics.EstimationSettlingTime_s = estimationSettlingTime;
metrics.AngleRMSLast1s_deg = sqrt(mean( ...
    rad2deg(x(1,lastSecondState)).^2));
metrics.ControlRMSLast1s_Nm = sqrt(mean( ...
    uActuator(lastSecondInput).^2));
metrics.FiniteDifferenceOmegaRMSE_rad_s = sqrt(mean( ...
    (x(2,postTransientState) ...
     - omegaFiniteDifference(postTransientState)).^2));
metrics.ControllerPoleRadius = max(abs(design.controllerPoles));
metrics.ObserverPoleRadius = max(abs(design.observerPoles));
metrics.ModelMassRatio = config.model.m/plant.m;
metrics.ModelLengthRatio = config.model.l/plant.l;
end

function entryTime = firstPermanentEntry(t, condition)
if isempty(condition) || ~condition(end)
    entryTime = Inf;
    return;
end
lastOutsideIndex = find(~condition, 1, 'last');
if isempty(lastOutsideIndex)
    entryTime = t(1);
elseif lastOutsideIndex == numel(condition)
    entryTime = Inf;
else
    entryTime = t(lastOutsideIndex + 1);
end
end

function metricsTable = resultsToTable(results)
rows = repmat(metricRow(results(1)), 1, numel(results));
for resultIndex = 1:numel(results)
    rows(resultIndex) = metricRow(results(resultIndex));
end
metricsTable = struct2table(rows);
end

function row = metricRow(result)
config = result.config;
metrics = result.metrics;
row.Case = string(config.name);
row.Group = string(config.group);
row.Mode = string(config.mode);
row.Ts_s = config.Ts;
row.Qtheta = config.Q(1,1);
row.Qomega = config.Q(2,2);
row.R = config.R;
row.Theta0_deg = rad2deg(config.x0(1));
row.Omega0_rad_s = config.x0(2);
row.TorqueLimit_Nm = max(abs([config.umin, config.umax]));
row.ModelMassRatio = metrics.ModelMassRatio;
row.ModelLengthRatio = metrics.ModelLengthRatio;
row.InputNoiseStd_Nm = config.inputDisturbanceStd;
row.MeasurementNoiseStd_deg = rad2deg(config.measurementNoiseStd);
row.SettlingTime_s = metrics.SettlingTime_s;
row.Settled = metrics.Settled;
row.FinalAngle_deg = metrics.FinalAngle_deg;
row.FinalOmega_rad_s = metrics.FinalOmega_rad_s;
row.MaxAngle_deg = metrics.MaxAngle_deg;
row.MaxTorque_Nm = metrics.MaxTorque_Nm;
row.ControlEffort_Nm2s = metrics.ControlEffort_Nm2s;
row.MaxPlantTorque_Nm = metrics.MaxPlantTorque_Nm;
row.PlantInputEffort_Nm2s = metrics.PlantInputEffort_Nm2s;
row.Saturation_pct = metrics.Saturation_pct;
row.AngleEstimateRMSE_deg = metrics.AngleEstimateRMSE_deg;
row.OmegaEstimateRMSE_rad_s = metrics.OmegaEstimateRMSE_rad_s;
row.EstimationSettlingTime_s = metrics.EstimationSettlingTime_s;
row.AngleRMSLast1s_deg = metrics.AngleRMSLast1s_deg;
row.ControlRMSLast1s_Nm = metrics.ControlRMSLast1s_Nm;
row.ControllerPoleRadius = metrics.ControllerPoleRadius;
row.ObserverPoleRadius = metrics.ObserverPoleRadius;
row.FiniteDifferenceOmegaRMSE_rad_s = ...
    metrics.FiniteDifferenceOmegaRMSE_rad_s;
end

function design = designDiscreteSystem(model, gravity, Ts, Q, R, C, ...
        requestedObserverPoles)
[Ac, Bc] = continuousLinearModel(model, gravity);
zohMatrix = expm([Ac, Bc; zeros(1,3)]*Ts);
Ad = zohMatrix(1:2, 1:2);
Bd = zohMatrix(1:2, 3);

assert(rank([Bd, Ad*Bd]) == 2, ...
    'The discrete model is not controllable.');
assert(rank([C; C*Ad]) == 2, ...
    'The discrete model is not observable from angle alone.');

[K, controllerPoles] = designDiscreteLQR(Ad, Bd, Q, R);
L = designObserverGain(Ad, C, requestedObserverPoles);
observerPoles = eig(Ad - L*C);

assert(all(abs(controllerPoles) < 1), ...
    'The discrete controller poles are not inside the unit circle.');
assert(all(abs(observerPoles) < 1), ...
    'The discrete observer poles are not inside the unit circle.');

design.Ac = Ac;
design.Bc = Bc;
design.Ad = Ad;
design.Bd = Bd;
design.C = C;
design.K = K;
design.L = L;
design.controllerPoles = controllerPoles;
design.observerPoles = observerPoles;
end

function [Ac, Bc] = continuousLinearModel(model, gravity)
Ac = [0, 1; ...
      model.m*gravity*model.l/model.J, -model.b/model.J];
Bc = [0; 1/model.J];
end

function dx = pendulumDynamics(~, x, u, parameters, gravity)
dx = [x(2); ...
      (parameters.m*gravity*parameters.l*sin(x(1)) ...
       - parameters.b*x(2) + u)/parameters.J];
end

function [K, closedLoopPoles] = designDiscreteLQR(A, B, Q, R)
if exist('dlqr', 'file') == 2
    [K, ~, closedLoopPoles] = dlqr(A, B, Q, R);
    return;
end

P = Q;
maximumIterations = 10000;
converged = false;
for iteration = 1:maximumIterations
    Kcurrent = (R + B'*P*B) \ (B'*P*A);
    Pnext = A'*P*A - A'*P*B*Kcurrent + Q;
    tolerance = 1e-12*max(1, norm(Pnext, 'fro'));
    if norm(Pnext - P, 'fro') <= tolerance
        P = Pnext;
        converged = true;
        break;
    end
    P = Pnext;
end
assert(converged, 'The discrete Riccati iteration did not converge.');
K = (R + B'*P*B) \ (B'*P*A);
closedLoopPoles = eig(A - B*K);
end

function L = designObserverGain(A, C, desiredPoles)
if exist('place', 'file') == 2
    L = place(A', C', desiredPoles)';
    return;
end

dualA = A';
dualB = C';
dualControllability = [dualB, dualA*dualB];
assert(rank(dualControllability) == 2, ...
    'The dual system is not controllable.');
desiredPolynomial = poly(desiredPoles);
polynomialOfA = dualA^2 + desiredPolynomial(2)*dualA ...
    + desiredPolynomial(3)*eye(2);
dualGain = ([0, 1]/dualControllability)*polynomialOfA;
L = dualGain';
end

function [time, data] = extractLoggedSignal(simulationOutput, names)
time = [];
data = [];
for nameIndex = 1:numel(names)
    candidate = [];
    try
        candidate = simulationOutput.get(names{nameIndex});
    catch
        % The variable is not a direct SimulationOutput property.
    end
    [time, data] = normalizeLoggedValue(candidate);
    if ~isempty(data)
        return;
    end
end

try
    logsout = simulationOutput.get('logsout');
    for nameIndex = 1:numel(names)
        try
            element = logsout.getElement(names{nameIndex});
            [time, data] = normalizeLoggedValue(element.Values);
            if ~isempty(data)
                return;
            end
        catch
            % Try the next possible signal name.
        end
    end
catch
    % No logsout dataset was available.
end
end

function [time, data] = normalizeLoggedValue(value)
time = [];
data = [];
if isempty(value)
    return;
end
if isa(value, 'Simulink.SimulationData.Signal')
    value = value.Values;
end
if isa(value, 'timeseries')
    time = double(value.Time(:));
    data = squeeze(double(value.Data));
    data = data(:);
elseif isstruct(value) && isfield(value, 'time') ...
        && isfield(value, 'signals')
    time = double(value.time(:));
    data = squeeze(double(value.signals.values));
    data = data(:);
end
end

function effort = loggedSignalEffort(time, data)
if numel(time) < 2 || numel(time) ~= numel(data)
    effort = NaN;
else
    effort = trapz(time(:), data(:).^2);
end
end

function saveSimulinkFigure(resultsDir, thetaTime, thetaData, ...
        thetaHatTime, thetaHatData, omegaTime, omegaData, ...
        omegaHatTime, omegaHatData, actuatorTorqueTime, ...
        actuatorTorqueData, plantTorqueTime, plantTorqueData)
fig = figure(5); clf(fig);
set(fig, 'Name', '5 - Simulink validation', 'NumberTitle', 'off', ...
    'Visible', 'on', 'Color', 'w', 'Position', [100, 100, 900, 700]);
layout = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(layout);
plot(ax, thetaTime, rad2deg(thetaData), 'LineWidth', 1.3, ...
    'DisplayName', 'True');
hold(ax, 'on');
if ~isempty(thetaHatData)
    plot(ax, thetaHatTime, rad2deg(thetaHatData), '--', ...
        'LineWidth', 1.1, 'DisplayName', 'Estimated');
end
grid(ax, 'on'); ylabel(ax, 'Angle [deg]');
title(ax, 'Final Simulink/Simscape validation');
legend(ax, 'Location', 'best');

ax = nexttile(layout);
if isempty(omegaData)
    text(ax, 0.5, 0.5, 'Angular velocity was not logged', ...
        'HorizontalAlignment', 'center');
    axis(ax, 'off');
else
    plot(ax, omegaTime, omegaData, 'LineWidth', 1.3, ...
        'DisplayName', 'True');
    hold(ax, 'on');
    if ~isempty(omegaHatData)
        plot(ax, omegaHatTime, omegaHatData, '--', ...
            'LineWidth', 1.1, 'DisplayName', 'Estimated');
    end
    grid(ax, 'on'); ylabel(ax, 'Velocity [rad/s]');
    legend(ax, 'Location', 'best');
end

ax = nexttile(layout);
if isempty(actuatorTorqueData) && isempty(plantTorqueData)
    text(ax, 0.5, 0.5, 'Torque was not logged', ...
        'HorizontalAlignment', 'center');
    axis(ax, 'off');
else
    hold(ax, 'on');
    if ~isempty(actuatorTorqueData)
        stairs(ax, actuatorTorqueTime, actuatorTorqueData, ...
            'LineWidth', 1.2, 'DisplayName', 'Actuator torque');
    end
    if ~isempty(plantTorqueData)
        stairs(ax, plantTorqueTime, plantTorqueData, '--', ...
            'LineWidth', 1.1, 'DisplayName', 'Total plant torque');
    end
    grid(ax, 'on'); ylabel(ax, 'Torque [N*m]'); xlabel(ax, 'Time [s]');
    legend(ax, 'Location', 'best');
end

exportgraphics(fig, fullfile(resultsDir, ...
    '05_simulink_validation.png'), 'Resolution', 160);
end
