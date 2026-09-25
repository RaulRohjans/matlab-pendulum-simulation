%% CONTROL AND INTELLIGENT SYSTEMS
% Project Part 1 - Stabilizing an Actuated Inverted Pendulum
%
% This script:
%   1.  Simulates the nonlinear pendulum (u = 0) with ode45        [Sec 1.1]
%   2.  Derives and analyses the linearized model                  [Sec 1.2]
%   3.  Discretizes the linear model with a zero-order hold        [Sec 1.3]
%   4.  Designs a discrete-time LQR                                [Sec 2]
%   5.  Evaluates it on the NONLINEAR pendulum:                    [Sec 2.2]
%         - linear vs nonlinear, with/without actuator saturation
%         - different Q and R
%         - different initial conditions
%         - input disturbance
%         - model uncertainty (controller model m.* ~= plant p.*)
%   6.  Designs a Luenberger observer (angle-only measurement)     [Sec 3]
%   7.  Observer experiments:                                      [Sec 3.2]
%         - wrong initial estimate, observer poles, measurement noise,
%           sampling time, finite-difference velocity estimate
%   8.  Runs the supplied Simulink/Simscape model
%
% Convention: theta = 0 is the UPRIGHT position, CCW positive.
% "Plant" simulations in MATLAB integrate the nonlinear ODE with ode45
% between samples while the input is held constant (ZOH), i.e. they
% mimic the sampled-data loop of the Simulink model.

clear; clc; close all;

% Choose which experiments to run (each one is a few seconds)
RUN.openloop    = true;
RUN.lqr         = true;
RUN.QR_sweep    = true;
RUN.IC_sweep    = true;
RUN.disturbance = true;
RUN.uncertainty = true;
RUN.observer    = true;
RUN.Ts_sweep    = true;
RUN.simulink    = true;


%% ================================================================
%  1. PHYSICAL PARAMETERS  (the "real" pendulum, used by Simscape)
%  ================================================================

p.m = 1.0;           % mass [kg]
p.l = 0.5;           % distance pivot -> centre of mass [m]
p.b = 0.05;          % viscous damping coefficient [N*m/(rad/s)]
g   = 9.81;          % gravitational acceleration [m/s^2]

p.theta0 = 35;       % initial angle [deg]
p.omega0 = 0;        % initial angular velocity. NOTE: the Revolute Joint in
                     % the Simscape model reads this value in deg/s.


%% ================================================================
%  2. MODEL PARAMETERS  (used for control/observer design)
%  ================================================================
% Identical to the physical parameters for now. Section 9 changes them
% on purpose to study model uncertainty.

m.m = p.m;
m.l = p.l;
m.b = p.b;
m.J = m.m * m.l^2;   % moment of inertia [kg*m^2]


%% ================================================================
%  3. SAMPLING TIME AND ACTUATOR LIMITS
%  ================================================================

Ts = 0.02;           % sampling time [s]

% Holding the pendulum at angle theta needs a torque of m*g*l*sin(theta)
% (2.8 N*m at 35 deg, 4.9 N*m at 90 deg), so the limit must exceed this.
umin = -10;          % [N*m]
umax =  10;          % [N*m]

Tsim   = 5;          % length of the closed-loop simulations [s]
tolSet = deg2rad(1); % settling band |theta| <= 1 deg

x_initial = [deg2rad(p.theta0); p.omega0];


%% ================================================================
%  4. PART 1.1 - NONLINEAR MODEL, OPEN LOOP (u = 0)
%  ================================================================
% J*theta_ddot = m*g*l*sin(theta) - b*theta_dot + u
% x1 = theta, x2 = theta_dot:
%     x1_dot = x2
%     x2_dot = (m*g*l*sin(x1) - b*x2 + u)/J

if RUN.openloop
    tspan   = [0 10];
    odeOpts = odeset('RelTol',1e-9,'AbsTol',1e-11);  % tight: the equilibrium
                                                     % is unstable, so solver
                                                     % errors grow quickly
    % [theta0 (rad), omega0 (rad/s)]
    ICs = [ 0              0
            1e-6           0
            deg2rad(1)     0
           -deg2rad(1)     0
            deg2rad(5)     0
            0              0.01 ];
    icNames = {'[0, 0]', '\theta_0 = 10^{-6} rad', '\theta_0 = +1°', ...
               '\theta_0 = -1°', '\theta_0 = 5°', '\omega_0 = 0.01 rad/s'};

    figure('Name','1.1 Open loop');
    for i = 1:size(ICs,1)
        [t_ol, x_ol] = ode45(@(t,x) nonlinear_pendulum(t,x,0,p,g), ...
                             tspan, ICs(i,:)', odeOpts);
        subplot(2,1,1); plot(t_ol, rad2deg(x_ol(:,1)), 'LineWidth',1.3); hold on;
        subplot(2,1,2); plot(rad2deg(x_ol(:,1)), x_ol(:,2), 'LineWidth',1.3); hold on;
    end
    subplot(2,1,1); grid on; yline(180,'k:'); yline(-180,'k:');
    xlabel('Time [s]'); ylabel('\theta [deg]');
    title('Nonlinear pendulum without control (\theta = 0 upright)');
    legend(icNames, 'Location','eastoutside');
    subplot(2,1,2); grid on;
    xlabel('\theta [deg]'); ylabel('\omega [rad/s]'); title('Phase portrait');
    legend(icNames, 'Location','eastoutside');
    % [0,0] is an equilibrium and stays there. Any perturbation grows, the
    % pendulum falls to the side of the perturbation and, due to damping,
    % settles at the hanging equilibrium theta = +-180 deg.
end


%% ================================================================
%  5. PART 1.2 - LINEARIZED MODEL AND STABILITY
%  ================================================================
% Around theta = 0, theta_dot = 0, u = 0 with sin(theta) ~ theta:
%     theta_ddot = (m*g*l*theta - b*theta_dot + u)/J

Ac = [ 0,               1
       m.m*g*m.l/m.J,  -m.b/m.J ];
Bc = [ 0
       1/m.J ];

continuous_poles = eig(Ac);

fprintf('\n=========== Continuous-time linear model ===========\n');
disp('Ac ='); disp(Ac);
disp('Bc ='); disp(Bc);
disp('Eigenvalues of Ac:'); disp(continuous_poles);
if any(real(continuous_poles) > 0)
    fprintf(['One eigenvalue is positive (%.2f 1/s): the upright equilibrium ' ...
             'of the linearization is UNSTABLE (saddle).\n'], max(real(continuous_poles)));
    fprintf(['By Lyapunov''s indirect method the upright equilibrium of the ' ...
             'nonlinear pendulum is therefore also unstable,\nconsistent with ' ...
             'Section 1.1. Small deviations grow roughly like exp(%.2f t).\n'], ...
             max(real(continuous_poles)));
end

if RUN.openloop
    % Linear vs nonlinear for a small initial angle
    x0s   = [deg2rad(1); 0];
    t_cmp = linspace(0, 2, 400);
    [~, x_nl] = ode45(@(t,x) nonlinear_pendulum(t,x,0,p,g), t_cmp, x0s, odeOpts);
    x_lin = zeros(2, numel(t_cmp));
    for k = 1:numel(t_cmp)
        x_lin(:,k) = expm(Ac*t_cmp(k)) * x0s;
    end
    figure('Name','1.2 Linear vs nonlinear');
    plot(t_cmp, rad2deg(x_nl(:,1)), 'LineWidth',1.5); hold on; grid on;
    plot(t_cmp, rad2deg(x_lin(1,:)), '--', 'LineWidth',1.5);
    ylim([-5 200]); xlabel('Time [s]'); ylabel('\theta [deg]');
    title('\theta_0 = 1°: the linear model is accurate only while \theta is small');
    legend('Nonlinear (ode45)', 'Linearized', 'Location','northwest');
end


%% ================================================================
%  6. PART 1.3 - DISCRETE-TIME MODEL (ZOH)
%  ================================================================
% The input is held constant between samples: x[k+1] = Ad x[k] + Bd u[k],
% with Ad = expm(Ac*Ts) and Bd = int_0^Ts expm(Ac*tau) d tau * Bc.
% In Simulink, the Zero-Order Hold blocks sample theta and theta_dot every
% Ts (measurements) and hold the controller output constant between
% samples (actuation), which is exactly what this model assumes.

n = size(Ac,1);
sys = expm([Ac Bc; zeros(size(Bc,2), n+size(Bc,2))] * Ts);
Ad = sys(1:n, 1:n);
Bd = sys(1:n, n+1:end);

fprintf('\n=========== Discrete-time model (Ts = %.3f s) ===========\n', Ts);
disp('Ad ='); disp(Ad);
disp('Bd ='); disp(Bd);
disp('Discrete open-loop poles eig(Ad):'); disp(eig(Ad));
disp('Check: exp(eig(Ac)*Ts) gives the same poles:'); disp(exp(continuous_poles*Ts));
% A pole with |z| > 1 means the discrete model is unstable, as expected.


%% ================================================================
%  7. PART 2 - DISCRETE-TIME LQR
%  ================================================================
% Minimize J = sum( x'Qx + u'Ru ).
% Q penalizes state deviation, R penalizes control effort.

Q = diag([10, 1]);
R = 0.1;

[K, ~, cl_poles] = dlqr(Ad, Bd, Q, R);

fprintf('\n=========== Discrete-time LQR ===========\n');
disp('Q ='); disp(Q);
disp('R ='); disp(R);
disp('K ='); disp(K);
disp('Closed-loop poles eig(Ad - Bd*K):'); disp(cl_poles);
disp('|closed-loop poles| (must be < 1):'); disp(abs(cl_poles));

% Nominal design, used by all experiments below and by Simulink
design = struct('Ad',Ad, 'Bd',Bd, 'K',K, 'C',[1 0], 'L',[0;0], 'Ts',Ts);


%% ================================================================
%  8. PART 2.2 - LQR EXPERIMENTS (full state measured)
%  ================================================================
% Base options for simulate_loop (see local functions at the end)
base = struct('umin',-Inf, 'umax',Inf, 'useObserver',false, ...
              'xhat0',[0;0], 'distAmp',0, 'noiseAmp',0, ...
              'linearPlant',false, 'seed',1);

% ---- 8a. Linear vs nonlinear plant, with/without saturation ----------
if RUN.lqr
    o = base;                  o.linearPlant = true;
    r_lin    = simulate_loop(p, g, design, x_initial, Tsim, o);
    o = base;
    r_nl     = simulate_loop(p, g, design, x_initial, Tsim, o);
    o = base;  o.umin = umin;  o.umax = umax;
    r_nl_sat = simulate_loop(p, g, design, x_initial, Tsim, o);

    figure('Name','2.2 LQR: linear vs nonlinear, saturation');
    subplot(2,1,1);
    plot(r_lin.t, rad2deg(r_lin.x(1,:)), 'LineWidth',1.5); hold on; grid on;
    plot(r_nl.t,  rad2deg(r_nl.x(1,:)),  'LineWidth',1.5);
    plot(r_nl_sat.t, rad2deg(r_nl_sat.x(1,:)), '--', 'LineWidth',1.5);
    ylabel('\theta [deg]');
    title(sprintf('LQR from \\theta_0 = %g°', p.theta0));
    legend('Linear plant', 'Nonlinear plant', 'Nonlinear + saturation');
    subplot(2,1,2);
    stairs(r_lin.t(1:end-1), r_lin.u, 'LineWidth',1.5); hold on; grid on;
    stairs(r_nl.t(1:end-1),  r_nl.u,  'LineWidth',1.5);
    stairs(r_nl_sat.t(1:end-1), r_nl_sat.u, '--', 'LineWidth',1.5);
    yline(umin,'k:'); yline(umax,'k:');
    xlabel('Time [s]'); ylabel('u [N m]');
    legend('Linear plant', 'Nonlinear plant', 'Nonlinear + saturation', ...
           'u_{min}', 'u_{max}');

    fprintf('\n=========== LQR on the nonlinear pendulum ===========\n');
    print_table({'linear plant','nonlinear','nonlinear + saturation'}, ...
        [performance(r_lin,tolSet), performance(r_nl,tolSet), ...
         performance(r_nl_sat,tolSet)]);
end

% ---- 8b. Effect of Q and R ------------------------------------------
if RUN.QR_sweep
    QR_cases = { diag([1 1]),     1
                 diag([10 1]),    0.1
                 diag([100 1]),   0.1
                 diag([10 1]),    10
                 diag([1000 10]), 1 };
    names = cell(1, size(QR_cases,1));
    res   = [];
    figure('Name','2.2 Q/R trade-off');
    for i = 1:size(QR_cases,1)
        Qi = QR_cases{i,1};  Ri = QR_cases{i,2};
        di = design;  di.K = dlqr(Ad, Bd, Qi, Ri);
        o = base;  o.umin = umin;  o.umax = umax;
        r = simulate_loop(p, g, di, x_initial, Tsim, o);
        res = [res, performance(r, tolSet)]; %#ok<AGROW>
        names{i} = sprintf('Q=diag(%g,%g), R=%g', Qi(1,1), Qi(2,2), Ri);
        subplot(2,1,1); plot(r.t, rad2deg(r.x(1,:)), 'LineWidth',1.3); hold on;
        subplot(2,1,2); stairs(r.t(1:end-1), r.u, 'LineWidth',1.3); hold on;
    end
    subplot(2,1,1); grid on; ylabel('\theta [deg]');
    title('Effect of Q and R (nonlinear plant, saturated)'); legend(names);
    subplot(2,1,2); grid on; ylabel('u [N m]'); xlabel('Time [s]');
    yline(umin,'k:'); yline(umax,'k:');
    fprintf('\n=========== Q/R sweep ===========\n');
    print_table(names, res);
    % Larger Q (or smaller R): faster settling and smaller angle, but larger
    % torques that hit the saturation. Larger R: gentler torques, slower
    % response and a larger excursion.
end

% ---- 8c. Different initial conditions (with and without limits) -----
% Gravity torque is at most m*g*l = 4.9 N*m, so with |u| <= 10 N*m the
% pendulum can always be pushed back up. A tight limit below m*g*l
% (umax_tight) shows that stabilization is then only possible from small
% enough initial angles.
if RUN.IC_sweep
    umax_tight  = 3;                           % [N*m], < m*g*l
    theta0_list = [10 35 45 60 90 120 150 175];   % deg
    names = arrayfun(@(a) sprintf('theta0 = %g deg', a), theta0_list, ...
                     'UniformOutput', false);
    res_free = [];  res_sat = [];  res_tight = [];
    figure('Name','2.2 Initial conditions');
    for i = 1:numel(theta0_list)
        x0i = [deg2rad(theta0_list(i)); 0];
        r1 = simulate_loop(p, g, design, x0i, Tsim, base);
        o = base;  o.umin = umin;  o.umax = umax;
        r2 = simulate_loop(p, g, design, x0i, Tsim, o);
        o = base;  o.umin = -umax_tight;  o.umax = umax_tight;
        r3 = simulate_loop(p, g, design, x0i, Tsim, o);
        res_free  = [res_free,  performance(r1,tolSet)]; %#ok<AGROW>
        res_sat   = [res_sat,   performance(r2,tolSet)]; %#ok<AGROW>
        res_tight = [res_tight, performance(r3,tolSet)]; %#ok<AGROW>
        subplot(1,3,1); plot(r1.t, rad2deg(r1.x(1,:)), 'LineWidth',1.2); hold on;
        subplot(1,3,2); plot(r2.t, rad2deg(r2.x(1,:)), 'LineWidth',1.2); hold on;
        subplot(1,3,3); plot(r3.t, rad2deg(r3.x(1,:)), 'LineWidth',1.2); hold on;
    end
    subplot(1,3,1); grid on; title('No actuator limits');
    xlabel('Time [s]'); ylabel('\theta [deg]');
    subplot(1,3,2); grid on; title(sprintf('|u| \\leq %g N m', umax));
    xlabel('Time [s]');
    subplot(1,3,3); grid on; title(sprintf('|u| \\leq %g N m', umax_tight));
    xlabel('Time [s]'); legend(names, 'Location','best');
    fprintf('\n=========== Initial conditions, no limits ===========\n');
    print_table(names, res_free);
    fprintf('\n=========== Initial conditions, |u| <= %g ===========\n', umax);
    print_table(names, res_sat);
    fprintf('\n=========== Initial conditions, |u| <= %g ===========\n', umax_tight);
    print_table(names, res_tight);
    % Without limits the LQR stabilizes even large angles, at the cost of
    % huge torques. Saturation makes the response slower. When the limit is
    % below m*g*l the pendulum can no longer be held beyond
    % asin(umax/(m*g*l)) and larger initial angles cannot be recovered -
    % the region of attraction is limited.
end

% ---- 8d. Disturbance on the control input ---------------------------
% u_plant[k] = sat(-K x[k]) + v[k], v uniform in [-0.05, 0.05] N*m
% (same as the Uniform Random Number block in the Simulink model).
if RUN.disturbance
    o = base;  o.umin = umin;  o.umax = umax;
    r0 = simulate_loop(p, g, design, x_initial, Tsim, o);
    o.distAmp = 0.05;
    r1 = simulate_loop(p, g, design, x_initial, Tsim, o);
    o.distAmp = 0.5;
    r2 = simulate_loop(p, g, design, x_initial, Tsim, o);

    figure('Name','2.2 Input disturbance');
    subplot(2,1,1);
    plot(r0.t, rad2deg(r0.x(1,:)), r1.t, rad2deg(r1.x(1,:)), ...
         r2.t, rad2deg(r2.x(1,:)), 'LineWidth',1.3); grid on;
    ylabel('\theta [deg]'); title('LQR with input disturbance');
    legend('no disturbance', '|v| \leq 0.05 N m', '|v| \leq 0.5 N m');
    subplot(2,1,2);   % zoom on the steady state
    plot(r1.t, rad2deg(r1.x(1,:)), r2.t, rad2deg(r2.x(1,:)), 'LineWidth',1.3);
    grid on; xlim([2 Tsim]); xlabel('Time [s]'); ylabel('\theta [deg]');
    title('Zoom: steady-state fluctuation');
    fprintf('\n=========== Input disturbance ===========\n');
    print_table({'none','|v|<=0.05','|v|<=0.5'}, ...
        [performance(r0,tolSet), performance(r1,tolSet), performance(r2,tolSet)]);
end


%% ================================================================
%  9. PART 2.2 - MODEL UNCERTAINTY
%  ================================================================
% The controller is designed with a WRONG model (mdl) and applied to the
% true plant p. (To do this in Simulink, change m.* in Section 2.)
if RUN.uncertainty
    unc = { 'nominal',          1.0, 1.0, 1.0
            'mass +30%',        1.3, 1.0, 1.0
            'mass -30%',        0.7, 1.0, 1.0
            'length +30%',      1.0, 1.3, 1.0
            'length -30%',      1.0, 0.7, 1.0
            'damping x5',       1.0, 1.0, 5.0
            'm -50%, l -30%',   0.5, 0.7, 1.0
            'length -70%',      1.0, 0.3, 1.0 };
    names = unc(:,1)';
    res = [];
    figure('Name','2.2 Model uncertainty');
    for i = 1:size(unc,1)
        mdl   = struct('m', unc{i,2}*p.m, 'l', unc{i,3}*p.l, 'b', unc{i,4}*p.b);
        d_unc = design_controller(mdl, g, Ts, Q, R, []);
        o = base;  o.umin = umin;  o.umax = umax;
        r = simulate_loop(p, g, d_unc, x_initial, Tsim, o);   % true plant p
        res = [res, performance(r, tolSet)]; %#ok<AGROW>
        plot(r.t, rad2deg(r.x(1,:)), 'LineWidth',1.3); hold on;
    end
    grid on; xlabel('Time [s]'); ylabel('\theta [deg]');
    title('Controller designed with a wrong model, applied to the true plant');
    legend(names);
    fprintf('\n=========== Model uncertainty ===========\n');
    print_table(names, res);
end


%% ================================================================
%  10. PART 3 - OBSERVER DESIGN
%  ================================================================
% Only theta is measured: y[k] = C x[k], C = [1 0].
% Observer (prediction form, matches the Simulink Unit Delay structure):
%     xhat[k+1] = Ad xhat[k] + Bd u[k] + L (y[k] - C xhat[k])
% Estimation error: e[k+1] = (Ad - L C) e[k]  -> choose eig(Ad - L C).
%
% Why the Unit Delay in Simulink: the observer computes xhat[k+1] from
% xhat[k]. The Unit Delay stores that output for one sample and feeds it
% back as xhat[k] at the next step (initialised with xhat0). Without it
% the observer output would depend on itself instantaneously (algebraic
% loop) - it is the discrete-time "memory" (state) of the observer.

C = [1 0];

fprintf('\n=========== Observer ===========\n');
fprintf('Rank of observability matrix: %d (need 2)\n', rank(obsv(Ad, C)));

% Observer poles: faster than the controller poles (|z| smaller), but not
% so fast that measurement noise is amplified. Rule of thumb: 3-10x faster.
observer_poles = [0.5 0.6];
L = place(Ad', C', observer_poles)';

disp('L ='); disp(L);
disp('Observer poles eig(Ad - L*C):');  disp(eig(Ad - L*C));
disp('Controller poles eig(Ad - Bd*K):'); disp(eig(Ad - Bd*K));
fprintf('Equivalent continuous-time observer poles: %.1f, %.1f 1/s\n', ...
        log(observer_poles)/Ts);
fprintf('Equivalent continuous-time controller poles: %.1f, %.1f 1/s\n', ...
        real(log(eig(Ad - Bd*K))/Ts));

design.C = C;
design.L = L;
obs_s = log(observer_poles)/Ts;   % same observer speed for other Ts

% Initial estimate for Simulink (Unit Delay initial condition)
xhat0 = [0; 0];     % deliberately wrong: true initial angle is p.theta0


%% ================================================================
%  11. PART 3.2 - OBSERVER EXPERIMENTS
%  ================================================================
if RUN.observer
    % ---- 11a. Output feedback with correct vs wrong initial estimate ----
    o = base;  o.umin = umin;  o.umax = umax;  o.useObserver = true;
    o.xhat0 = x_initial;           r_ok  = simulate_loop(p, g, design, x_initial, Tsim, o);
    o.xhat0 = [0; 0];              r_bad = simulate_loop(p, g, design, x_initial, Tsim, o);
    o.xhat0 = [0; 0]; o.useObserver = false;   % full state, observer in parallel
    r_fs = simulate_loop(p, g, design, x_initial, Tsim, o);

    figure('Name','3.2 Observer: initial estimate');
    subplot(3,1,1);
    plot(r_bad.t, rad2deg(r_bad.x(1,:)), 'LineWidth',1.5); hold on; grid on;
    plot(r_bad.t, rad2deg(r_bad.xhat(1,:)), '--', 'LineWidth',1.5);
    plot(r_fs.t,  rad2deg(r_fs.x(1,:)), ':', 'LineWidth',1.5);
    ylabel('\theta [deg]'); legend('actual', 'estimate', 'full-state LQR');
    title('Observer-based LQR, wrong initial estimate xhat0 = [0; 0]');
    subplot(3,1,2);
    plot(r_bad.t, r_bad.x(2,:), 'LineWidth',1.5); hold on; grid on;
    plot(r_bad.t, r_bad.xhat(2,:), '--', 'LineWidth',1.5);
    ylabel('\omega [rad/s]'); legend('actual', 'estimate');
    subplot(3,1,3);
    stairs(r_bad.t(1:end-1), r_bad.u, 'LineWidth',1.5); hold on; grid on;
    stairs(r_ok.t(1:end-1),  r_ok.u,  '--', 'LineWidth',1.5);
    xlabel('Time [s]'); ylabel('u [N m]');
    legend('wrong xhat0', 'correct xhat0');

    fprintf('\n=========== Observer-based LQR ===========\n');
    print_table({'full state','observer, correct xhat0','observer, wrong xhat0'}, ...
        [performance(r_fs,tolSet), performance(r_ok,tolSet), performance(r_bad,tolSet)]);

    % ---- 11b. Observer pole choice, with and without measurement noise ---
    pole_sets = { [0.98 0.99], [0.8 0.85], [0.5 0.6], [0.1 0.15] };
    noise     = deg2rad(0.5);     % uniform measurement noise amplitude
    names = cellfun(@(z) sprintf('z = [%g %g]', z), pole_sets, 'UniformOutput', false);
    res_clean = [];  res_noisy = [];
    rmse = zeros(numel(pole_sets), 2);
    figure('Name','3.2 Observer poles');
    for i = 1:numel(pole_sets)
        di = design;  di.L = place(Ad', C', pole_sets{i})';
        o = base;  o.umin = umin;  o.umax = umax;  o.useObserver = true;
        r1 = simulate_loop(p, g, di, x_initial, Tsim, o);
        o.noiseAmp = noise;
        r2 = simulate_loop(p, g, di, x_initial, Tsim, o);
        res_clean = [res_clean, performance(r1,tolSet)]; %#ok<AGROW>
        res_noisy = [res_noisy, performance(r2,tolSet)]; %#ok<AGROW>
        idx = find(r2.t >= 2);       % after the transient
        rmse(i,:) = [rmsv(r2.x(2,idx) - r2.xhat(2,idx)), ...
                     rmsv(r1.x(2,:)   - r1.xhat(2,:))];
        subplot(2,1,1); plot(r1.t, r1.x(2,:) - r1.xhat(2,:), 'LineWidth',1.2); hold on;
        subplot(2,1,2); plot(r2.t, r2.x(2,:) - r2.xhat(2,:), 'LineWidth',1.2); hold on;
    end
    subplot(2,1,1); grid on; ylabel('\omega - \omega_{hat} [rad/s]');
    title('Velocity estimation error, no measurement noise'); legend(names);
    subplot(2,1,2); grid on; ylabel('\omega - \omega_{hat} [rad/s]');
    xlabel('Time [s]'); ylim([-1 1]);
    title(sprintf('With measurement noise |w| \\leq %.1f°', rad2deg(noise)));
    fprintf('\n=========== Observer poles, no noise ===========\n');
    print_table(names, res_clean);
    fprintf('\n=========== Observer poles, measurement noise ===========\n');
    print_table(names, res_noisy);
    disp('RMS velocity estimation error [rad/s]:');
    disp(array2table(rmse, 'RowNames', names, ...
         'VariableNames', {'noisy_after_2s','noise_free_total'}));
    % Very slow observer (slower than the controller): the wrong initial
    % estimate takes seconds to die out and the loop does not settle.
    % Very fast observer: quick convergence but it passes measurement noise
    % straight into the velocity estimate and the torque (higher effort).

    % ---- 11c. Alternative: finite-difference velocity estimate ----------
    o = base;  o.umin = umin;  o.umax = umax;  o.useObserver = true;
    o.noiseAmp = noise;
    r = simulate_loop(p, g, design, x_initial, Tsim, o);
    w_fd = [0, diff(r.y)/Ts];                  % (y[k]-y[k-1])/Ts
    tk   = r.t(1:end-1);
    figure('Name','3.2 Observer vs finite difference');
    plot(r.t, r.x(2,:), 'k', 'LineWidth',1.5); hold on; grid on;
    stairs(r.t, r.xhat(2,:), 'LineWidth',1.2);
    stairs(tk, w_fd, 'LineWidth',1.0);
    xlabel('Time [s]'); ylabel('\omega [rad/s]');
    legend('actual', 'observer', 'finite difference (y_k - y_{k-1})/T_s');
    title('Velocity estimation from a noisy angle measurement');
    idx  = find(tk >= 2);
    fprintf('RMS error after 2 s: observer %.3f rad/s, finite difference %.3f rad/s\n', ...
        rmsv(r.x(2,idx) - r.xhat(2,idx)), rmsv(r.x(2,idx) - w_fd(idx)));
    % Finite differencing amplifies noise by ~1/Ts and ignores the model.
    % The observer uses the model (and u) to filter the measurement.
    % Other options: low-pass filtered derivative, Kalman filter.
end


%% ================================================================
%  12. PART 3.2 - SAMPLING TIME
%  ================================================================
% Controller and observer are redesigned for each Ts with the same Q, R
% and the same continuous-time observer poles obs_s.
if RUN.Ts_sweep
    Ts_list = [0.005 0.02 0.05 0.1 0.2 0.3];
    names = arrayfun(@(T) sprintf('Ts = %g s', T), Ts_list, 'UniformOutput', false);
    res_fs = [];  res_obs = [];
    figure('Name','3.2 Sampling time');
    for i = 1:numel(Ts_list)
        di = design_controller(m, g, Ts_list(i), Q, R, obs_s);
        o = base;  o.umin = umin;  o.umax = umax;
        r1 = simulate_loop(p, g, di, x_initial, Tsim, o);         % full state
        o.useObserver = true;
        r2 = simulate_loop(p, g, di, x_initial, Tsim, o);         % observer
        res_fs  = [res_fs,  performance(r1,tolSet)]; %#ok<AGROW>
        res_obs = [res_obs, performance(r2,tolSet)]; %#ok<AGROW>
        subplot(1,2,1); plot(r1.t, rad2deg(r1.x(1,:)), 'LineWidth',1.2); hold on;
        subplot(1,2,2); plot(r2.t, rad2deg(r2.x(1,:)), 'LineWidth',1.2); hold on;
    end
    subplot(1,2,1); grid on; title('Full-state LQR');
    xlabel('Time [s]'); ylabel('\theta [deg]');
    subplot(1,2,2); grid on; title('Observer-based LQR (xhat0 = 0)');
    xlabel('Time [s]'); legend(names);
    fprintf('\n=========== Sampling time, full state ===========\n');
    print_table(names, res_fs);
    fprintf('\n=========== Sampling time, observer ===========\n');
    print_table(names, res_obs);
    % Larger Ts: the plant drifts further between samples (unstable pole
    % exp(4.33*Ts) grows), the controller reacts later and the observer has
    % fewer measurements. The observer-based loop degrades first (bigger
    % overshoot, more effort) and fails around Ts = 0.3 s.
end


%% ================================================================
%  13. RUN THE SUPPLIED SIMULINK MODEL
%  ================================================================
% Variables used by actuated_pendulum.slx:
%   p.m, p.l, p.b, p.theta0, p.omega0, Ts, K, umin, umax, Ad, Bd, L, C, xhat0
% Check the wiring before running:
%   - controller output -> Zero-Order Hold -> pendulum "input torque"
%   - observer inputs: u = controller output, y = sampled theta
%   - controller input = xhat (Unit Delay) for output feedback,
%     or the sampled [theta; theta_dot] (Mux) for full-state LQR
%   - disturbance: Add block between the controller and the pendulum,
%     with the Uniform Random Number block as the second input
% Controller block:  u = min(max(-K*x, umin), umax);
% Observer block:    xhatnext = Ad*xhat + Bd*u + L*(y - C*xhat);

% Make sure the nominal design is in the workspace
Ad = design.Ad;  Bd = design.Bd;  K = design.K;  L = design.L;  C = design.C;

if RUN.simulink
    if exist('actuated_pendulum', 'file') == 4
        out = sim('actuated_pendulum');

        t_s  = out.theta_sim.Time;   th_s = squeeze(out.theta_sim.Data);
        w_s  = squeeze(out.omega_sim.Data);
        tu_s = out.u_sim.Time;       u_s  = squeeze(out.u_sim.Data);
        th_h = squeeze(out.theta_hat.Data);   t_h = out.theta_hat.Time;
        w_h  = squeeze(out.omega_hat.Data);

        figure('Name','Simulink / Simscape result');
        subplot(3,1,1);
        plot(t_s, rad2deg(th_s), 'LineWidth',1.5); hold on; grid on;
        stairs(t_h, rad2deg(th_h), '--', 'LineWidth',1.2);
        ylabel('\theta [deg]'); legend('measured', 'estimate');
        title('Simscape pendulum with discrete LQR + observer');
        subplot(3,1,2);
        plot(t_s, w_s, 'LineWidth',1.5); hold on; grid on;
        stairs(t_h, w_h, '--', 'LineWidth',1.2);
        ylabel('\omega [rad/s]'); legend('actual', 'estimate');
        subplot(3,1,3);
        stairs(tu_s, u_s, 'LineWidth',1.5); hold on; grid on;
        yline(umin,'k:'); yline(umax,'k:');
        xlabel('Time [s]'); ylabel('u [N m]');
    else
        warning('actuated_pendulum.slx not found on the MATLAB path - skipped.');
    end
end


%% ================================================================
%  14. SAVE IMPORTANT PARAMETERS
%  ================================================================
save('pendulum_parameters.mat', 'p', 'm', 'g', 'Ts', 'umin', 'umax', ...
     'Ac', 'Bc', 'Ad', 'Bd', 'Q', 'R', 'K', 'C', 'L', 'xhat0');


%% ================================================================
%  LOCAL FUNCTIONS
%  ================================================================

function dx = nonlinear_pendulum(~, x, u, par, g)
% Nonlinear pendulum dynamics, x = [theta; theta_dot], theta = 0 upright.
% par has fields m, l, b (J is computed here).
    J  = par.m * par.l^2;
    dx = [ x(2)
           (par.m*g*par.l*sin(x(1)) - par.b*x(2) + u) / J ];
end

function d = design_controller(mdl, g, Ts, Q, R, obs_s)
% Linearize, discretize (ZOH) and design LQR (+ observer if obs_s given)
% from the model parameters mdl (fields m, l, b).
    J  = mdl.m * mdl.l^2;
    Ac = [0 1; mdl.m*g*mdl.l/J, -mdl.b/J];
    Bc = [0; 1/J];
    sd = c2d(ss(Ac, Bc, eye(2), zeros(2,1)), Ts, 'zoh');
    d.Ad = sd.A;  d.Bd = sd.B;  d.Ts = Ts;  d.C = [1 0];
    d.K  = dlqr(d.Ad, d.Bd, Q, R);
    if isempty(obs_s)
        d.L = [0; 0];
    else
        d.L = place(d.Ad', d.C', exp(obs_s*Ts))';
    end
end

function r = simulate_loop(plant, g, d, x0, T, o)
% Sampled-data closed loop: controller/observer at Ts, input held (ZOH),
% plant integrated with ode45 (or the discrete linear model).
%   o.useObserver : controller uses xhat (true) or the true state (false);
%                   the observer always runs, so it can be compared.
%   o.distAmp     : uniform input disturbance amplitude [N*m]
%   o.noiseAmp    : uniform measurement noise amplitude [rad]
    Ts = d.Ts;  N = round(T/Ts);
    x  = zeros(2, N+1);  xh = zeros(2, N+1);
    u  = zeros(1, N);    y  = zeros(1, N);
    x(:,1) = x0;  xh(:,1) = o.xhat0;

    rng(o.seed);
    v = o.distAmp  * (2*rand(1,N) - 1);
    w = o.noiseAmp * (2*rand(1,N) - 1);
    odeOpts = odeset('RelTol',1e-6, 'AbsTol',1e-8);

    for k = 1:N
        y(k) = d.C*x(:,k) + w(k);                     % measurement

        if o.useObserver, xc = xh(:,k); else, xc = x(:,k); end
        u(k) = min(max(-d.K*xc, o.umin), o.umax);      % LQR + saturation
        u_plant = u(k) + v(k);                         % input disturbance

        if o.linearPlant
            x(:,k+1) = d.Ad*x(:,k) + d.Bd*u_plant;
        else
            [~, xs]  = ode45(@(t,s) nonlinear_pendulum(t, s, u_plant, plant, g), ...
                             [0 Ts], x(:,k), odeOpts);
            x(:,k+1) = xs(end,:)';
        end

        xh(:,k+1) = d.Ad*xh(:,k) + d.Bd*u(k) + d.L*(y(k) - d.C*xh(:,k));
    end
    r = struct('t', (0:N)*Ts, 'x', x, 'xhat', xh, 'u', u, 'y', y, 'Ts', Ts);
end

function M = performance(r, tol)
% Settling time = first time after which |theta| STAYS within tol.
    th  = r.x(1,:);
    out = find(abs(th) > tol, 1, 'last');
    if isempty(out)
        ts = 0;
    elseif out == numel(th)
        ts = NaN;                                      % never settled
    else
        ts = r.t(out + 1);
    end
    M.stabilized     = ~isnan(ts);
    M.settling_time  = ts;
    M.max_angle_deg  = max(abs(rad2deg(th)));
    M.max_u          = max(abs(r.u));
    M.control_effort = sum(r.u.^2) * r.Ts;             % approx. int u^2 dt
end

function v = rmsv(e)
    v = sqrt(mean(e.^2));
end

function print_table(names, M)
    Tbl = struct2table(M(:));
    Tbl.Properties.RowNames = names(:);
    disp(Tbl);
end