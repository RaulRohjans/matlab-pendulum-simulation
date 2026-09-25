%% Physical parameters
%%%%% These are the parameters used by the Simulink model and they
%%%%% correspond to the properties of the real physical system

p.m = 1.0;        % kg
p.l = 0.5;        % m, distance pivot -> point mass
p.b = 0.05;       % N*m/(rad/s)
g = 9.81;       % m/s^2

p.theta0 = 35;  %degrees
p.omega0 = 0;   %rad/s

%% Model parameters
%%%%%%% These are the parameters you can use for your model. Currently,
%%%%%%% they are set equal to the values of the physical system but you can
%%%%%%% change and observe how your controller and observer behave under model
%%%%%%% uncertainty

m.m = p.m;        % kg
m.l = p.l;        % m, distance pivot -> point mass
m.b = p.b;       % N*m/(rad/s)
m.J = m.m*m.l^2;


%% Sampling
%%%%%% Controller works with measurements obtained at regular intervals.
%%%%%% The control inputs are in discrete time as well. They are then
%%%%%% applied to the system in zero-order-hold manner.
%%%%%% Sampling period: every Ts seconds
Ts = 0.02;

%% Actuator limits
% Chosen so the default initial condition (theta0 = 35 deg) can be
% recovered comfortably (peak unsaturated |u| ~ 24 N*m for this design),
% while still being tight enough that lowering these later (Section 2.2,
% "Actuator limits") visibly changes the closed-loop behavior.
umin = -30;
umax =  30;

%% Linear model
%%%%%% Here you should set the matrices Ac and Bc which correspond to
%%%%%% the linearized system
%
% From J*thetaddot = m*g*l*theta - b*thetadot + u (small-angle linearization
% of Eq. (1) about theta=0, thetadot=0, u=0), with x = [theta; thetadot]:

Ac = [0                 1;
      g/m.l           -m.b/m.J];
Bc = [0; 1/m.J];

%% Discrete-time model
%%%%% Disretize the linear apprroximation; you can use MATLAB functions
%%%%% from control systtem toolbox such as c2d or compute matrix
%%%%% exponentials

sysc = ss(Ac, Bc, eye(2), 0);
sysd = c2d(sysc, Ts, 'zoh');
Ad = sysd.A;
Bd = sysd.B;

%% LQR
%%%%% Design LQR

% Weights: penalize angle error much more than angular velocity, and
% keep R small enough to allow the (large) unsaturated torque needed to
% recover from theta0 = 35 deg. See the report for the Q/R trade-off
% discussion (settling time vs. control effort vs. saturation).
Q = diag([100, 1]);
R = 0.05;

K = dlqr(Ad, Bd, Q, R);

% Sanity check: closed-loop poles must satisfy abs(pole) < 1
disp('Closed-loop poles (LQR):'); disp(eig(Ad - Bd*K));

%% Observer
%%%% Design an observer using pole placement
C = [1 0];

% Deliberately WRONG initial estimate (true initial angle is 35 deg =
% p.theta0, but the observer starts assuming the pendulum is upright) --
% this is what lets you see estimator convergence in the simulation.
xhat0 = [0; 0];

% Observer error poles chosen well inside the unit circle and faster
% (smaller magnitude) than the LQR closed-loop poles above, so the state
% estimate settles before it meaningfully affects control performance.
obs_poles = [0.3, 0.35];
L = place(Ad', C', obs_poles)';

disp('Observer error poles:'); disp(eig(Ad - L*C));

%% Simulate the actuated pendulum
%%% The following command runs the simulation.
%%% You should see an animation of the pendulum and out is a struct that
%%% contains the outputs from the simulation (for example, measured and
%%% estimated values, control input, or any other thing you add)

out = sim("actuated_pendulum");

%%%%%%
%%% Add relevant plots etc.

figure('Name', 'Closed-loop response');
subplot(3,1,1);
plot(theta_sim.Time, rad2deg(theta_sim.Data), 'b', 'DisplayName', '\theta (true)'); hold on;
plot(theta_hat.Time, rad2deg(theta_hat.Data), 'r--', 'DisplayName', '\theta hat (estimate)');
ylabel('\theta [deg]'); legend show; grid on;
title('Angle: true vs. estimated');

subplot(3,1,2);
plot(omega_sim.Time, omega_sim.Data, 'b', 'DisplayName', '\thetadot (true)'); hold on;
plot(omega_hat.Time, omega_hat.Data, 'r--', 'DisplayName', '\thetadot hat (estimate)');
ylabel('\thetadot [rad/s]'); legend show; grid on;

subplot(3,1,3);
plot(u_sim.Time, u_sim.Data); hold on;
yline(umax, 'k--'); yline(umin, 'k--');
ylabel('u [N m]'); xlabel('t [s]'); grid on;
title('Control input (dashed lines = actuator limits)');
