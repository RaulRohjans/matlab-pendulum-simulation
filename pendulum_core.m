%% Actuated inverted pendulum - core implementation
% Control and Intelligent Systems, Project Part 1
%
% This script contains only the components needed before the experiments:
%   1. nonlinear continuous-time pendulum model (ode45),
%   2. linearization about the upright equilibrium,
%   3. exact zero-order-hold discretization,
%   4. discrete-time LQR controller,
%   5. angle-only discrete observer,
%   6. one nonlinear closed-loop reference simulation and result plots,
%   7. one run of the working supplied Simulink/Simscape model.
%
% The variables p, m, Ts, K, umin, umax, Ad, Bd, C, L, and xhat0 remain
% in the workspace so that the Simulink model can use them later.
%
% Control System Toolbox is used when available. Small two-state fallbacks
% are included for dlqr/place so the equation-based checks run in base MATLAB.
% The final supplied-model run needs Simulink and Simscape Multibody.

clear; clc; close all;

%% 1. Physical and model parameters
% p: physical plant parameters used by the nonlinear simulation.
p.m = 1.0;          % mass [kg]
p.l = 0.5;          % pivot-to-centre-of-mass distance [m]
p.b = 0.05;         % viscous damping [N*m*s/rad]
p.J = p.m*p.l^2;    % moment of inertia [kg*m^2]
p.theta0 = 35;      % initial angle [deg]
p.omega0 = 0;       % initial angular velocity [rad/s]

g = 9.81;           % gravitational acceleration [m/s^2]

% m: parameters used to design the controller and observer. Keep these
% equal to p until the model-uncertainty experiment is performed.
m.m = p.m;
m.l = p.l;
m.b = p.b;
m.J = m.m*m.l^2;

Ts = 0.02;          % controller sampling period [s]
Tend = 5;           % simulation duration [s]
umin = -10;         % minimum actuator torque [N*m]
umax =  10;         % maximum actuator torque [N*m]

runSimulinkModel = true; % set false to run only the equation-based checks

%% 2. Validate the nonlinear physical model with u = 0
% State x = [theta; omega], where theta = 0 is the upright position.
odeOptions = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
freeDynamics = @(t, x) pendulumDynamics(t, x, 0, p, g);

[tEquilibrium, xEquilibrium] = ode45( ...
    freeDynamics, [0 Tend], [0; 0], odeOptions);
[tPositive, xPositive] = ode45( ...
    freeDynamics, [0 Tend], [deg2rad(1); 0], odeOptions);
[tNegative, xNegative] = ode45( ...
    freeDynamics, [0 Tend], [-deg2rad(1); 0], odeOptions);

figure('Name', 'Open-loop nonlinear pendulum');
subplot(2,1,1);
plot(tEquilibrium, rad2deg(xEquilibrium(:,1)), 'LineWidth', 1.4);
hold on;
plot(tPositive, rad2deg(xPositive(:,1)), 'LineWidth', 1.4);
plot(tNegative, rad2deg(xNegative(:,1)), 'LineWidth', 1.4);
yline(180, 'k:', 'HandleVisibility', 'off');
yline(-180, 'k:', 'HandleVisibility', 'off');
grid on;
xlabel('Time [s]');
ylabel('Angle from upright [deg]');
title('Nonlinear pendulum without control (u = 0)');
legend('Exact equilibrium', '+1 deg perturbation', ...
    '-1 deg perturbation', 'Location', 'best');

% The phase portrait makes the two symmetric falls from the unstable
% upright equilibrium visible. With damping, both trajectories spiral
% toward a hanging equilibrium at +180 or -180 degrees.
subplot(2,1,2);
plot(rad2deg(xPositive(:,1)), xPositive(:,2), 'LineWidth', 1.4);
hold on;
plot(rad2deg(xNegative(:,1)), xNegative(:,2), 'LineWidth', 1.4);
grid on;
xlabel('Angle from upright [deg]');
ylabel('Angular velocity [rad/s]');
title('Phase portrait');
legend('+1 deg perturbation', '-1 deg perturbation', 'Location', 'best');

%% 3. Linear continuous-time model about theta = 0
% J*theta_ddot = m*g*l*theta - b*theta_dot + u
Ac = [0,                    1;
      m.m*g*m.l/m.J, -m.b/m.J];
Bc = [0; 1/m.J];

openLoopPoles = eig(Ac);
fprintf('Continuous-time open-loop poles:\n');
disp(openLoopPoles);
fprintf('The positive pole confirms that the upright equilibrium is unstable.\n\n');

%% 4. Exact zero-order-hold discrete model
% expm([Ac Bc; 0 0]*Ts) = [Ad Bd; 0 1]
zohMatrix = expm([Ac, Bc; zeros(1,3)]*Ts);
Ad = zohMatrix(1:2, 1:2);
Bd = zohMatrix(1:2, 3);

%% 5. Discrete-time LQR controller
Q = diag([100, 1]); % prioritize angle error over angular velocity
R = 1;              % penalize control effort
[K, controllerPoles] = designDiscreteLQR(Ad, Bd, Q, R);

controllabilityMatrix = [Bd, Ad*Bd];
assert(rank(controllabilityMatrix) == 2, ...
    'The discrete model is not controllable.');
assert(all(abs(controllerPoles) < 1), ...
    'The designed discrete-time controller is not stable.');

%% 6. Angle-only observer
% Measurement: y[k] = C*x[k]. The predictor observer is
% xhat[k+1] = Ad*xhat[k] + Bd*u[k] + L*(y[k] - C*xhat[k]).
% In Simulink, the Unit Delay stores xhat[k] while the observer computes
% xhat[k+1], which gives the discrete loop memory and avoids an algebraic loop.
C = [1, 0];
xhat0 = [0; 0];
observerPoles = [0.55, 0.60]; % faster than the controller poles
L = designObserverGain(Ad, C, observerPoles);

observabilityMatrix = [C; C*Ad];
assert(rank(observabilityMatrix) == 2, ...
    'The discrete model is not observable from the angle measurement.');
assert(all(abs(eig(Ad - L*C)) < 1), ...
    'The designed discrete-time observer is not stable.');

fprintf('Ad =\n'); disp(Ad);
fprintf('Bd =\n'); disp(Bd);
fprintf('LQR gain K =\n'); disp(K);
fprintf('Controller poles =\n'); disp(controllerPoles);
fprintf('Observer gain L =\n'); disp(L);
fprintf('Observer poles =\n'); disp(eig(Ad - L*C));

%% 7. Nonlinear sampled-data closed-loop simulation
% At each sample the controller uses only xhat. The resulting torque is
% saturated and held constant while ode45 advances the nonlinear plant.
t = 0:Ts:Tend;
numberOfSamples = numel(t);

x = zeros(2, numberOfSamples);
xhat = zeros(2, numberOfSamples);
y = zeros(1, numberOfSamples);
u = zeros(1, numberOfSamples - 1);

x(:,1) = [deg2rad(p.theta0); p.omega0];
xhat(:,1) = xhat0;

for k = 1:numberOfSamples - 1
    y(k) = C*x(:,k);                    % angle is the only measurement
    uRequested = -K*xhat(:,k);          % discrete LQR control law
    u(k) = min(max(uRequested, umin), umax);

    innovation = y(k) - C*xhat(:,k);
    xhat(:,k+1) = Ad*xhat(:,k) + Bd*u(k) + L*innovation;

    heldInputDynamics = @(time, state) ...
        pendulumDynamics(time, state, u(k), p, g);
    [~, stateOverSample] = ode45( ...
        heldInputDynamics, [t(k), t(k+1)], x(:,k), odeOptions);
    x(:,k+1) = stateOverSample(end,:)';
end
y(end) = C*x(:,end);

assert(abs(x(1,end)) < deg2rad(0.5) && abs(x(2,end)) < 0.05, ...
    'The nominal nonlinear closed loop did not settle near the upright state.');
fprintf('Final angle: %.4f deg\n', rad2deg(x(1,end)));
fprintf('Final angular velocity: %.4f rad/s\n', x(2,end));
fprintf('Maximum applied torque: %.4f N*m\n', max(abs(u)));

%% 8. Plot the closed-loop response
figure('Name', 'Closed-loop response');

subplot(3,1,1);
plot(t, rad2deg(x(1,:)), 'LineWidth', 1.4);
hold on;
plot(t, rad2deg(xhat(1,:)), '--', 'LineWidth', 1.2);
grid on;
ylabel('Angle [deg]');
title('LQR control using only the measured angle');
legend('True', 'Estimated', 'Location', 'best');

subplot(3,1,2);
plot(t, x(2,:), 'LineWidth', 1.4);
hold on;
plot(t, xhat(2,:), '--', 'LineWidth', 1.2);
grid on;
ylabel('Angular velocity [rad/s]');
legend('True', 'Estimated', 'Location', 'best');

subplot(3,1,3);
stairs(t(1:end-1), u, 'LineWidth', 1.4, ...
    'DisplayName', 'Applied torque');
hold on;
yline(umax, 'k--', 'DisplayName', 'Limits');
yline(umin, 'k--', 'HandleVisibility', 'off');
grid on;
xlabel('Time [s]');
ylabel('Torque [N*m]');
legend('show', 'Location', 'best');

%% 9. Run the working supplied Simulink/Simscape model
% The Karthik attempt contains the correctly wired controller, observer,
% Zero-Order Hold, and Unit Delay. Load it by absolute path so this script
% works regardless of MATLAB's current folder.
projectRoot = fileparts(mfilename('fullpath'));
modelFile = fullfile(projectRoot, ...
    'karthik simple pendulm', 'actuated_pendulum.slx');

if runSimulinkModel
    if exist('sim', 'file') ~= 2
        warning(['Simulink is not installed, so the supplied-model run was ' ...
                 'skipped. The equation-based simulation completed.']);
    elseif ~isfile(modelFile)
        warning('The supplied model was not found: %s', modelFile);
    else
        [~, modelName] = fileparts(modelFile);
        load_system(modelFile);

        % The supplied geometry treats p.l as a full visual length and puts
        % its point mass at p.l/2. The PDF instead defines p.l as the
        % pivot-to-COM distance. Override only the simulation copy of p so
        % the Simscape COM is at the intended distance without editing the SLX.
        pForSimscape = p;
        pForSimscape.l = 2*p.l;
        simulationInput = Simulink.SimulationInput(modelName);
        simulationInput = simulationInput.setVariable('p', pForSimscape);
        simulationInput = simulationInput.setModelParameter( ...
            'StopTime', num2str(Tend));
        out = sim(simulationInput);

        % ReturnWorkspaceOutputs is enabled in the model, so logged signals
        % are properties of the SimulationOutput object out.
        try
            figure('Name', 'Simulink/Simscape closed-loop response');
            subplot(2,1,1);
            plot(out.theta_sim.Time, rad2deg(out.theta_sim.Data), ...
                'LineWidth', 1.4);
            hold on;
            plot(out.theta_hat.Time, rad2deg(out.theta_hat.Data), '--', ...
                'LineWidth', 1.2);
            grid on;
            ylabel('Angle [deg]');
            title('Supplied nonlinear Simscape model');
            legend('True', 'Estimated', 'Location', 'best');

            subplot(2,1,2);
            plot(out.omega_sim.Time, out.omega_sim.Data, 'LineWidth', 1.4);
            hold on;
            plot(out.omega_hat.Time, out.omega_hat.Data, '--', ...
                'LineWidth', 1.2);
            grid on;
            xlabel('Time [s]');
            ylabel('Angular velocity [rad/s]');
            legend('True', 'Estimated', 'Location', 'best');
        catch plotError
            warning('Simulation ran, but its logged signals could not be plotted: %s', ...
                plotError.message);
        end
    end
end

%% Local functions
function dx = pendulumDynamics(~, x, u, parameters, gravity)
% Nonlinear dynamics: J*theta_ddot = m*g*l*sin(theta) - b*theta_dot + u.
dx = [x(2);
      (parameters.m*gravity*parameters.l*sin(x(1)) ...
       - parameters.b*x(2) + u)/parameters.J];
end

function [K, closedLoopPoles] = designDiscreteLQR(A, B, Q, R)
% Use Control System Toolbox when present; otherwise solve the DARE by
% fixed-point iteration (sufficient for this small stabilizable system).
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

assert(converged, ...
    'The discrete Riccati iteration did not converge.');
K = (R + B'*P*B) \ (B'*P*A);
closedLoopPoles = eig(A - B*K);
end

function L = designObserverGain(A, C, desiredPoles)
% Use Control System Toolbox when present. The fallback applies Ackermann's
% formula to the dual two-state system (A', C').
if exist('place', 'file') == 2
    L = place(A', C', desiredPoles)';
    return;
end

dualA = A';
dualB = C';
dualControllability = [dualB, dualA*dualB];
assert(rank(dualControllability) == 2, ...
    'The dual system is not controllable, so observer poles cannot be placed.');

desiredPolynomial = poly(desiredPoles);
polynomialOfA = dualA^2 ...
    + desiredPolynomial(2)*dualA ...
    + desiredPolynomial(3)*eye(2);
dualGain = ([0, 1]/dualControllability)*polynomialOfA;
L = dualGain';
end
