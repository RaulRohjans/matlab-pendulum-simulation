%% pendulum_project.m
% =========================================================================
% Control and Intelligent Systems - Project Part 1
% Stabilizing an actuated inverted pendulum
%
% ONE FILE for the whole project. It runs every task and experiment from
% the assignment in plain MATLAB (no Simulink needed to run it):
%
%   Part 1  Model          1.1 ode45 simulation   1.2 linearization
%                          1.3 discretization (ZOH)
%   Part 2  LQR            2.1 design K
%                          2.2 initial conditions, Q/R, actuator limits,
%                              input disturbance, model uncertainty
%   Part 3  Observer       3.1 design L
%                          3.2 observer poles, wrong initial estimate,
%                              sampling time, measurement noise,
%                              alternative velocity estimators
%   Final   Demo: only the angle is measured + disturbances
%   Export  Sets K, umin, umax, L, ... in the Workspace for Simulink
%
% HOW TO USE
%   * Put this file in your project folder and press "Run" (F5), or run
%     one section at a time with "Run Section" (Ctrl+Enter).
%   * Sections marked [SETUP] must run before the experiments after them.
%   * Figures are saved as .png in the folder "figures".
%   * Needs: Control System Toolbox (c2d, dlqr, place, ctrb, obsv).
%   * Requires MATLAB R2016b or newer (local functions at the end).
%
% IMPORTANT: replace the parameter values in section 0 with the values
% from the .m file supplied with your Simulink model.
%
% HOW THE SIMULATION WORKS (same idea as the Simulink model)
%   Every Ts seconds:  measure y[k] -> estimate x_hat[k] -> u[k] = -K x_hat[k]
%   -> saturate -> add disturbance -> ZOH: hold u constant for Ts while
%   ode45 integrates the NONLINEAR pendulum (with the p. parameters).
% =========================================================================

clear; clc; close all;
rng(0);                         % fixed random seed -> same noise every run
if ~exist('figures', 'dir'), mkdir('figures'); end

%% 0. Parameters [SETUP] ==================================================
% p. = "real" pendulum (used to simulate the physical system)
p.m = 1.0;      % mass [kg]                        <-- replace
p.l = 0.5;      % axis -> center of mass [m]       <-- replace
p.b = 0.1;      % viscous damping [N m s/rad]      <-- replace
p.g = 9.81;     % gravity [m/s^2]

% m. = model used to DESIGN the controller and observer.
% Identical to p. for now; changed in the model uncertainty experiment.
m = p;

Ts   = 0.01;    % sampling period [s]
umin = -2;      % actuator torque limits [N m]     <-- your choice
umax =  2;

Tend = 5;       % simulation length for experiments [s]

% Default design weights (changed in the Q/R experiment)
Q = diag([10, 1]);   % weight on [theta, dtheta]
R = 1;               % weight on u
obs_speed = 4;       % observer poles = (controller poles)^obs_speed

%% 1.1 Physical model: nonlinear simulation with ode45, u = 0 ============
% State x = [x1; x2] = [theta; dtheta]
%   dx1/dt = x2
%   dx2/dt = (m g l sin(x1) - b x2 + u) / J,   J = m l^2
f0 = @(t, x) pendulum_ode(x, 0, p);           % u = 0

x0_list = [ 0       0;       % exactly at the equilibrium
            1e-3    0;       % small angle perturbations
           -1e-3    0;
            1e-2    0;
            0       1e-2 ];  % small velocity perturbation

figure('Name', '1.1 Nonlinear pendulum, u = 0'); hold on; grid on;
for i = 1:size(x0_list, 1)
    [t, y] = ode45(f0, [0 5], x0_list(i,:)');
    plot(t, y(:,1), 'LineWidth', 1.2, 'DisplayName', ...
        sprintf('\\theta_0 = %g, \\omega_0 = %g', x0_list(i,1), x0_list(i,2)));
end
yline(pi, 'k--', 'hanging down (\theta = \pi)', 'HandleVisibility', 'off');
yline(-pi, 'k--', 'HandleVisibility', 'off');
xlabel('time [s]'); ylabel('\theta [rad]'); legend('Location', 'best');
title('1.1  Nonlinear pendulum without control (u = 0)');
save_fig('1_1_nonlinear_u0');
% Expected: theta0 = 0 stays at 0 (it is an equilibrium). Any tiny
% perturbation makes the pendulum fall and swing around theta = +-pi.

%% 1.2 Linear approximation [SETUP] ======================================
% sin(theta) ~ theta near theta = 0:
%   J ddtheta = m g l theta - b dtheta + u
[Ac, Bc] = linear_model(m);

fprintf('\n===== 1.2 Linear model =====\n');
disp('Ac ='); disp(Ac);
disp('Bc ='); disp(Bc);
eig_c = eig(Ac);
fprintf('Eigenvalues of Ac: %s\n', mat2str(eig_c', 4));
fprintf('-> one eigenvalue is POSITIVE, so the upright equilibrium is UNSTABLE.\n');
fprintf('   A small error grows like exp(%.2f t): doubles every %.3f s.\n', ...
        max(eig_c), log(2)/max(eig_c));

% Compare linear vs nonlinear for a small and a larger initial angle
figure('Name', '1.2 Linear vs nonlinear');
th0s = [0.01, 0.3];
for i = 1:2
    subplot(1, 2, i); hold on; grid on;
    [t, y]  = ode45(f0, [0 1.5], [th0s(i); 0]);
    [tl, yl] = ode45(@(t, x) Ac*x, [0 1.5], [th0s(i); 0]);
    plot(t, y(:,1), 'LineWidth', 1.2); plot(tl, yl(:,1), '--', 'LineWidth', 1.2);
    ylim([-0.5 4]); xlabel('time [s]'); ylabel('\theta [rad]');
    legend('nonlinear', 'linear', 'Location', 'northwest');
    title(sprintf('\\theta_0 = %.2f rad', th0s(i)));
end
sgtitle('1.2  The linear model is accurate only near \theta = 0');
save_fig('1_2_linear_vs_nonlinear');

%% 1.3 Discrete-time model with zero-order hold [SETUP] ==================
[Ad, Bd] = discretize(Ac, Bc, Ts);

fprintf('\n===== 1.3 Discrete model (Ts = %g s) =====\n', Ts);
disp('Ad ='); disp(Ad);
disp('Bd ='); disp(Bd);
fprintf('|eig(Ad)| = %s  (> 1 means unstable)\n', mat2str(abs(eig(Ad))', 4));
fprintf('Check: eig(Ad) = exp(eig(Ac)*Ts)? max difference = %.1e\n', ...
        max(abs(sort(eig(Ad)) - sort(exp(eig_c*Ts)))));

% Same result "by hand" with the matrix exponential (augmented-matrix trick):
%   expm([Ac Bc; 0 0]*Ts) = [Ad Bd; 0 1]
M = expm([Ac Bc; zeros(1, 3)] * Ts);
fprintf('Check with expm: max difference = %.1e\n', ...
        max(max(abs(M(1:2, :) - [Ad Bd]))));
fprintf('Controllable? rank(ctrb(Ad,Bd)) = %d (need 2)\n', rank(ctrb(Ad, Bd)));

%% 2.1 Discrete-time LQR design [SETUP] ==================================
% Minimize  J = sum( x' Q x + u' R u ),  control law u[k] = -K x[k]
[K, ~, cl_poles] = dlqr(Ad, Bd, Q, R);

fprintf('\n===== 2.1 LQR =====\n');
fprintf('K = %s\n', mat2str(K, 4));
fprintf('Closed-loop poles: %s\n', mat2str(cl_poles.', 4));
fprintf('|poles| = %s  (all < 1 -> stable)\n', mat2str(abs(cl_poles)', 4));

% Base settings for all simulations (each experiment changes a few fields)
base.p = p;  base.Ts = Ts;  base.Tend = Tend;
base.K = K;  base.umin = umin;  base.umax = umax;
base.use_sat   = false;         % actuator limits off until 2.2c
base.estimator = 'full';        % 'full' | 'observer' | 'fd' | 'filtered'
base.L = [];  base.Ad = Ad;  base.Bd = Bd;
base.xhat0 = [0; 0];
base.dist_std  = 0;             % input disturbance  (std, N m)
base.noise_std = 0;             % measurement noise  (std, rad)
base.x0 = [0.2; 0];

r = simulate_loop(base);
plot_run(r, '2.1  LQR, full state, nonlinear pendulum, \theta_0 = 0.2 rad');
save_fig('2_1_lqr_basic');
print_metrics('2.1 basic run', r);

%% 2.2a Evaluation: different initial states =============================
fprintf('\n===== 2.2a Initial conditions (no actuator limits) =====\n');
th0_list = [0.1 0.3 0.5 0.8 1.0 1.2 1.5 2.0 2.5 3.0];
figure('Name', '2.2a Initial conditions'); hold on; grid on;
fprintf('%8s %8s %10s %10s\n', 'theta0', 'stable', 'max|u|', 'settle[s]');
for th0 = th0_list
    c = base;  c.x0 = [th0; 0];
    r = simulate_loop(c);
    fprintf('%8.2f %8d %10.2f %10.2f\n', th0, r.success, r.max_u, r.t_settle);
    plot(r.t, r.x(1,:), 'DisplayName', sprintf('\\theta_0 = %.1f', th0));
end
xlabel('time [s]'); ylabel('\theta [rad]'); legend('Location', 'eastoutside');
title('2.2a  LQR from different initial angles (no torque limit)');
save_fig('2_2a_initial_conditions');
% Without torque limits, LQR recovers even from large angles in simulation,
% but it needs very large torques (see max|u|). The design is based on the
% linear model, so far from theta = 0 there is no stability guarantee.

%% 2.2b Choosing Q and R ==================================================
fprintf('\n===== 2.2b Q and R experiments (theta0 = 0.3 rad) =====\n');
QR_list = { diag([1 1]),    1,    'Q=diag(1,1),   R=1';
            diag([10 1]),   1,    'Q=diag(10,1),  R=1   (default)';
            diag([100 1]),  1,    'Q=diag(100,1), R=1';
            diag([10 1]),   0.1,  'Q=diag(10,1),  R=0.1';
            diag([10 1]),   10,   'Q=diag(10,1),  R=10';
            diag([10 10]),  1,    'Q=diag(10,10), R=1' };
fprintf('%-32s %10s %10s %10s %10s\n', 'setting', 'settle[s]', ...
        'max|th|', 'effort', 'max|u|');
figure('Name', '2.2b Q and R');
for i = 1:size(QR_list, 1)
    Ki = dlqr(Ad, Bd, QR_list{i,1}, QR_list{i,2});
    c = base;  c.K = Ki;  c.x0 = [0.3; 0];
    r = simulate_loop(c);
    fprintf('%-32s %10.2f %10.3f %10.3f %10.2f\n', QR_list{i,3}, ...
            r.t_settle, r.max_theta, r.effort, r.max_u);
    subplot(2,1,1); hold on; plot(r.t, r.x(1,:), 'DisplayName', QR_list{i,3});
    subplot(2,1,2); hold on; stairs(r.tk, r.u, 'DisplayName', QR_list{i,3});
end
subplot(2,1,1); grid on; ylabel('\theta [rad]'); legend('Location', 'eastoutside');
title('2.2b  Effect of Q and R');
subplot(2,1,2); grid on; ylabel('u [N m]'); xlabel('time [s]');
legend('Location', 'eastoutside'); xlim([0 2]);
save_fig('2_2b_Q_R');
% Trade-off: larger Q (or smaller R) -> faster, smaller angle, MORE torque.
%            larger R (or smaller Q) -> gentler torque, SLOWER response.
% effort = sum(u^2)*Ts (energy-like measure of control use).

%% 2.2c Actuator limits (saturation) ======================================
fprintf('\n===== 2.2c Actuator limits [%g, %g] N m =====\n', umin, umax);
fprintf('Gravity torque m*g*l = %.2f N m. With |u| <= %.2f the motor can hold\n', ...
        p.m*p.g*p.l, umax);
fprintf('the pendulum only up to about asin(umax/(m g l)) = %.2f rad.\n', ...
        asin(min(1, umax/(p.m*p.g*p.l))));
fprintf('%8s %14s %14s\n', 'theta0', 'stable(no lim)', 'stable(limit)');
figure('Name', '2.2c Saturation');
th0_sat = [0.1 0.2 0.3 0.4 0.5 0.7];
for th0 = th0_sat
    c = base;  c.x0 = [th0; 0];
    r_free = simulate_loop(c);
    c.use_sat = true;
    r_sat = simulate_loop(c);
    fprintf('%8.2f %14d %14d\n', th0, r_free.success, r_sat.success);
    subplot(2,1,1); hold on; plot(r_sat.t, r_sat.x(1,:), ...
        'DisplayName', sprintf('\\theta_0 = %.1f', th0));
    subplot(2,1,2); hold on; stairs(r_sat.tk, r_sat.u, ...
        'DisplayName', sprintf('\\theta_0 = %.1f', th0));
end
subplot(2,1,1); grid on; ylabel('\theta [rad]'); ylim([-4 4]);
legend('Location', 'eastoutside'); title('2.2c  LQR with actuator limits');
subplot(2,1,2); grid on; ylabel('u [N m]'); xlabel('time [s]');
legend('Location', 'eastoutside');
save_fig('2_2c_saturation');
% With saturation the controller cannot give the torque LQR "wants".
% For large angles gravity wins and the pendulum falls: the region of
% initial states that can be stabilized becomes much smaller.

%% 2.2d Input disturbance: u[k] = -K x[k] + v[k] =========================
fprintf('\n===== 2.2d Input disturbance =====\n');
dist_list = [0 0.05 0.2 0.5];      % standard deviation of v [N m]
figure('Name', '2.2d Input disturbance'); hold on; grid on;
fprintf('%10s %12s %12s %8s\n', 'std(v)', 'RMS th(>2s)', 'max|th|', 'stable');
for d = dist_list
    c = base;  c.use_sat = true;  c.dist_std = d;  c.x0 = [0.1; 0];
    r = simulate_loop(c);
    fprintf('%10.2f %12.4f %12.3f %8d\n', d, r.rms_theta_end, r.max_theta, r.success);
    plot(r.t, r.x(1,:), 'DisplayName', sprintf('std(v) = %.2f N m', d));
end
xlabel('time [s]'); ylabel('\theta [rad]'); legend;
title('2.2d  LQR (with limits) under input disturbance');
save_fig('2_2d_input_disturbance');
% The controller keeps rejecting the disturbance, so theta stays small but
% keeps moving around 0. Bigger noise -> bigger deviations.

%% 2.2e Model uncertainty: design with m., simulate with p. ==============
fprintf('\n===== 2.2e Model uncertainty (real system p. is fixed) =====\n');
unc_list = { 'm.m', [0.5 0.8 1.2 1.5 2.0];
             'm.l', [0.5 0.8 1.2 1.5 2.0];
             'm.b', [0 0.5 2 5] };
fprintf('%-6s %8s %10s %10s %10s\n', 'param', 'factor', 'stable', ...
        'settle[s]', 'max|u|');
figure('Name', '2.2e Model uncertainty');
for j = 1:size(unc_list, 1)
    subplot(3, 1, j); hold on; grid on;
    for fac = unc_list{j, 2}
        mw = p;                                     % wrong model
        field = unc_list{j, 1}(3:end);
        mw.(field) = p.(field) * fac;
        [Acw, Bcw] = linear_model(mw);
        [Adw, Bdw] = discretize(Acw, Bcw, Ts);
        c = base;  c.K = dlqr(Adw, Bdw, Q, R);  c.use_sat = true;
        c.x0 = [0.2; 0];
        r = simulate_loop(c);
        fprintf('%-6s %8.2f %10d %10.2f %10.2f\n', unc_list{j,1}, fac, ...
                r.success, r.t_settle, r.max_u);
        plot(r.t, r.x(1,:), 'DisplayName', sprintf('%s x %.1f', unc_list{j,1}, fac));
    end
    ylabel('\theta [rad]'); legend('Location', 'eastoutside');
end
xlabel('time [s]'); sgtitle('2.2e  Controller designed with a wrong model');
save_fig('2_2e_model_uncertainty');
% LQR is quite robust: moderate model errors still stabilize the pendulum,
% but performance changes (slower / more oscillation). Very wrong models
% (e.g. mass or length much too small) can make it fail.

%% 3.1 Observer design [SETUP] ===========================================
% Only the angle is measured: y[k] = C x[k],  C = [1 0].
% Observer (predict with the model, then correct with the measurement):
%   x_pred[k] = Ad x_hat[k-1] + Bd u[k-1]          <- u[k-1] from Unit Delay
%   x_hat[k]  = x_pred[k] + L (y[k] - C x_pred[k])
% Estimation error: e[k] = (Ad - L C Ad) e[k-1]  -> choose L so that the
% eigenvalues of (Ad - L C Ad) are inside the unit circle and FAST.
C = [1 0];
fprintf('\n===== 3.1 Observer =====\n');
fprintf('Observable? rank(obsv(Ad,C)) = %d (need 2) -> dtheta CAN be estimated\n', ...
        rank(obsv(Ad, C)));

L = design_observer(Ad, C, cl_poles, obs_speed);
fprintf('Observer poles = (controller poles)^%d = %s\n', obs_speed, ...
        mat2str(eig(Ad - L*C*Ad).', 4));
fprintf('L = %s\n', mat2str(L', 4));

base.L = L;
c = base;  c.estimator = 'observer';  c.xhat0 = [0; 0];  c.x0 = [0.2; 0];
r = simulate_loop(c);
plot_run(r, '3.1  LQR + observer (only \theta measured)');
save_fig('3_1_observer_basic');
print_metrics('3.1 observer run', r);
% Why the Unit Delay? u[k] is computed from x_hat[k], and x_hat[k] needs an
% input. If it needed u[k] itself, the signal would depend on itself in the
% same time step (algebraic loop). The model step uses u[k-1], which is
% already known, so the Unit Delay provides exactly that.

%% 3.2a Observer poles and wrong initial estimate ========================
fprintf('\n===== 3.2a Observer pole speed (wrong initial estimate) =====\n');
speeds = [1 2 4 8];
fprintf('%8s %22s %14s %8s\n', 'speed', '|observer poles|', 'RMS est.err', 'stable');
figure('Name', '3.2a Observer poles');
for s = speeds
    Ls = design_observer(Ad, C, cl_poles, s);
    c = base;  c.estimator = 'observer';  c.L = Ls;
    c.x0 = [0.2; 0];  c.xhat0 = [0; 1];          % WRONG initial estimate
    r = simulate_loop(c);
    fprintf('%8d %22s %14.4f %8d\n', s, mat2str(abs(eig(Ad - Ls*C*Ad))', 3), ...
            r.rms_est_err, r.success);
    subplot(2,1,1); hold on;
    plot(r.tk, r.x_k(2,:) - r.xhat(2,:), 'DisplayName', sprintf('speed %d', s));
    subplot(2,1,2); hold on;
    plot(r.t, r.x(1,:), 'DisplayName', sprintf('speed %d', s));
end
subplot(2,1,1); grid on; xlim([0 1.5]); ylabel('\omega - \omega_{hat} [rad/s]');
legend; title('3.2a  Velocity estimation error (initial estimate is wrong)');
subplot(2,1,2); grid on; ylabel('\theta [rad]'); xlabel('time [s]'); legend;
save_fig('3_2a_observer_poles');
% Faster observer poles (closer to 0) -> the error disappears faster, but
% the first error peak can be bigger (large L reacts strongly at the start).
% Rule of thumb: observer 2-10x faster than the controller, so the
% controller gets a good estimate quickly. Too fast -> large L, which
% amplifies measurement noise (see 3.2c).

%% 3.2b Sampling time ====================================================
fprintf('\n===== 3.2b Sampling time =====\n');
Ts_list = [0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5];
fprintf('%8s %12s %12s %14s %14s\n', 'Ts', 'full: stable', 'obs: stable', ...
        'full: settle', 'obs: settle');
figure('Name', '3.2b Sampling time');
for i = 1:numel(Ts_list)
    Tsi = Ts_list(i);
    [Adi, Bdi] = discretize(Ac, Bc, Tsi);
    [Ki, ~, cpi] = dlqr(Adi, Bdi, Q, R);         % redesign for this Ts
    Li = design_observer(Adi, C, cpi, obs_speed);
    c = base;  c.Ts = Tsi;  c.K = Ki;  c.L = Li;  c.Ad = Adi;  c.Bd = Bdi;
    c.x0 = [0.2; 0];  c.use_sat = true;
    r_full = simulate_loop(c);                   % full state
    c.estimator = 'observer';  c.xhat0 = [0; 0];
    r_obs = simulate_loop(c);                    % observer
    fprintf('%8.3f %12d %12d %14.2f %14.2f\n', Tsi, r_full.success, ...
            r_obs.success, r_full.t_settle, r_obs.t_settle);
    subplot(2,1,1); hold on; plot(r_full.t, r_full.x(1,:), ...
        'DisplayName', sprintf('Ts = %g', Tsi));
    subplot(2,1,2); hold on; plot(r_obs.t, r_obs.x(1,:), ...
        'DisplayName', sprintf('Ts = %g', Tsi));
end
subplot(2,1,1); grid on; ylabel('\theta'); title('3.2b  Full state measured');
legend('Location', 'eastoutside'); ylim([-1 1]);
subplot(2,1,2); grid on; ylabel('\theta'); xlabel('time [s]');
title('Only \theta measured (observer)'); legend('Location', 'eastoutside');
ylim([-1 1]);
save_fig('3_2b_sampling_time');
% Larger Ts -> the controller reacts later and holds u constant for longer.
% Between samples the unstable pendulum drifts. At some Ts it fails.
% The observer version fails earlier (with the placeholder parameters: at
% Ts = 0.4 s with the observer, while full state still works at 0.5 s),
% because the estimate also gets worse with fewer samples.

% Full-state control with the observer running IN PARALLEL (only estimates)
c = base;  c.estimator = 'full';  c.run_observer_parallel = true;
c.x0 = [0.2; 0];  c.xhat0 = [0; 1];
r = simulate_loop(c);
figure('Name', '3.2b Observer in parallel');
subplot(2,1,1); plot(r.tk, r.x_k(1,:), r.tk, r.xhat(1,:), '--'); grid on;
ylabel('\theta'); legend('true', 'estimate');
title('Controller uses the TRUE state; observer only runs in parallel');
subplot(2,1,2); plot(r.tk, r.x_k(2,:), r.tk, r.xhat(2,:), '--'); grid on;
ylabel('\omega'); xlabel('time [s]'); legend('true', 'estimate');
save_fig('3_2b_observer_parallel');

%% 3.2c Measurement noise: y[k] = C x[k] + w[k] ==========================
fprintf('\n===== 3.2c Measurement noise =====\n');
noise_list = [0 0.001 0.005 0.01];               % std of w [rad]
speed_list = [2 4 8];
fprintf('%10s %8s %14s %12s %10s\n', 'std(w)', 'speed', 'RMS err(w)', ...
        'RMS u(>2s)', 'stable');
figure('Name', '3.2c Measurement noise');
k_plot = 0;
for s = speed_list
    Ls = design_observer(Ad, C, cl_poles, s);
    for w = noise_list
        c = base;  c.estimator = 'observer';  c.L = Ls;  c.use_sat = true;
        c.noise_std = w;  c.x0 = [0.1; 0];  c.xhat0 = c.x0;  % correct start
        r = simulate_loop(c);
        fprintf('%10.3f %8d %14.4f %12.4f %10d\n', w, s, r.rms_est_err, ...
                r.rms_u_end, r.success);
        if w == 0.005
            k_plot = k_plot + 1;
            subplot(numel(speed_list), 1, k_plot);
            stairs(r.tk, r.u); grid on; ylabel('u [N m]');
            title(sprintf('std(w) = %.3f rad, observer speed %d', w, s));
        end
    end
end
xlabel('time [s]'); sgtitle('3.2c  Control signal under measurement noise');
save_fig('3_2c_measurement_noise');
% The noise enters the estimate through L. A faster observer (bigger L)
% trusts the measurement more -> noisier velocity estimate -> noisier
% torque. Slower observer = smoother, but reacts later. Trade-off again.

%% 3.2d Alternative velocity estimators ==================================
% Compare: (1) observer, (2) finite difference (theta[k]-theta[k-1])/Ts,
%          (3) finite difference + first-order low-pass filter
fprintf('\n===== 3.2d Alternative estimators (with measurement noise) =====\n');
methods = {'observer', 'fd', 'filtered'};
names   = {'Observer', 'Finite difference', 'Filtered finite difference'};
figure('Name', '3.2d Alternative estimators');
fprintf('%-28s %16s %12s %8s\n', 'method', 'RMS vel. error', 'RMS u(>2s)', 'stable');
for i = 1:3
    c = base;  c.estimator = methods{i};  c.use_sat = true;
    c.noise_std = 0.002;  c.x0 = [0.1; 0];  c.xhat0 = c.x0;  % fair start
    c.filter_tau = 0.05;                          % low-pass time constant [s]
    r = simulate_loop(c);
    fprintf('%-28s %16.4f %12.4f %8d\n', names{i}, r.rms_est_err, ...
            r.rms_u_end, r.success);
    subplot(3,1,i);
    plot(r.tk, r.x_k(2,:), 'k', r.tk, r.xhat(2,:), 'r'); grid on;
    ylabel('\omega [rad/s]'); legend('true', 'estimate'); title(names{i});
end
xlabel('time [s]'); sgtitle('3.2d  Velocity estimates, std(w) = 0.002 rad');
save_fig('3_2d_alternative_estimators');
% Finite difference: simple, no model needed, but divides noise by Ts ->
% very noisy. The filter reduces noise but adds delay. The observer uses
% the model AND the known input, so it is usually the most accurate.
% (A Kalman filter = observer with L chosen optimally for given noise.)

%% FINAL DEMO: only the angle measured + input and measurement noise =====
fprintf('\n===== FINAL DEMO =====\n');
c = base;  c.estimator = 'observer';  c.use_sat = true;
c.dist_std = 0.05;  c.noise_std = 0.002;
c.x0 = [0.2; 0];  c.xhat0 = [0; 0];  c.Tend = 10;
r = simulate_loop(c);
plot_run(r, ['FINAL: LQR + observer, only \theta measured, ' ...
             'actuator limits, input + measurement noise']);
save_fig('final_demo');
print_metrics('final demo', r);

%% EXPORT: variables for the Simulink model ==============================
% Run this section before running the Simulink model. The MATLAB Function
% blocks read these values from the Workspace (set them as "Parameter"
% in the block's Symbols / Edit Data window).
m = p;                                    % back to the correct model
[Ac, Bc] = linear_model(m);
[Ad, Bd] = discretize(Ac, Bc, Ts);
[K, ~, cl_poles] = dlqr(Ad, Bd, Q, R);
L = design_observer(Ad, C, cl_poles, obs_speed);
xhat0 = [0; 0];
fprintf('\nWorkspace ready for Simulink: p, m, Ts, K, umin, umax, Ad, Bd, C, L, xhat0\n');
fprintf('Figures saved in folder: %s\n', fullfile(pwd, 'figures'));

% ---- Code to paste into the Simulink MATLAB Function blocks -------------
% Controller block:
%     function u = controller(x, K, umin, umax)
%         u = -K * x;
%         u = min(max(u, umin), umax);   % saturation (2.2c)
%     end
%
% Observer block (u_prev comes from the Unit Delay block):
%     function xhat = observer(y, u_prev, Ad, Bd, C, L, xhat0)
%         persistent xh
%         if isempty(xh), xh = xhat0; end
%         x_pred = Ad*xh + Bd*u_prev;          % predict
%         xh = x_pred + L*(y - C*x_pred);      % correct
%         xhat = xh;
%     end
% -------------------------------------------------------------------------


%% ========================================================================
%  LOCAL FUNCTIONS (MATLAB needs these at the END of a script file)
%  ========================================================================

function dx = pendulum_ode(x, u, p)
% Nonlinear pendulum:  J ddtheta = m g l sin(theta) - b dtheta + u
J  = p.m * p.l^2;
dx = [ x(2);
       (p.m*p.g*p.l*sin(x(1)) - p.b*x(2) + u) / J ];
end

function [Ac, Bc] = linear_model(m)
% Linearization around theta = 0, dtheta = 0, u = 0 (sin(theta) ~ theta)
J  = m.m * m.l^2;
Ac = [ 0,        1;
       m.g/m.l, -m.b/J ];
Bc = [ 0;
       1/J ];
end

function [Ad, Bd] = discretize(Ac, Bc, Ts)
% Zero-order hold discretization
sysd = c2d(ss(Ac, Bc, eye(2), 0), Ts, 'zoh');
Ad = sysd.A;  Bd = sysd.B;
end

function L = design_observer(Ad, C, cl_poles, speed)
% Current-estimator observer: error dynamics e[k] = (Ad - L*C*Ad) e[k-1].
% Observer poles = controller poles ^ speed  (z^n is n times "faster",
% because z = exp(s*Ts) -> z^n = exp(n*s*Ts)).
obs_poles = cl_poles(:) .^ speed;
if isreal(obs_poles) && abs(obs_poles(1) - obs_poles(2)) < 1e-6
    obs_poles(2) = obs_poles(2) * 0.99;   % place() needs distinct poles
end
L = place(Ad', (C*Ad)', obs_poles)';
end

function r = simulate_loop(c)
% Sampled-data closed loop: discrete controller + NONLINEAR pendulum.
% Between samples, u is held constant (ZOH) and ode45 integrates the plant.
if ~isfield(c, 'run_observer_parallel'), c.run_observer_parallel = false; end
if ~isfield(c, 'filter_tau'),            c.filter_tau = 0.05;            end

N  = round(c.Tend / c.Ts);
Cm = [1 0];
x  = c.x0(:);
xh = c.xhat0(:);                 % observer state
u_prev = 0;
y_prev = x(1);                   % for finite differences
w_filt = 0;                      % filtered velocity
a_filt = c.Ts / (c.filter_tau + c.Ts);   % low-pass coefficient

tk = (0:N-1) * c.Ts;
x_k = zeros(2, N);  xhat = zeros(2, N);  u_log = zeros(1, N);
t_cell = cell(1, N);  x_cell = cell(1, N);
opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
failed = false;

for k = 1:N
    % ---- measurement -------------------------------------------------
    y = Cm*x + c.noise_std*randn;

    % ---- state estimate ---------------------------------------------
    need_obs = strcmp(c.estimator, 'observer') || c.run_observer_parallel;
    if need_obs
        x_pred = c.Ad*xh + c.Bd*u_prev;          % predict (uses u[k-1])
        xh     = x_pred + c.L*(y - Cm*x_pred);   % correct with y[k]
    end
    switch c.estimator
        case 'full'
            x_used = [y; x(2)];              % angle (with noise) + true velocity
            if ~c.run_observer_parallel, xh = x_used; end
        case 'observer'
            x_used = xh;
        case 'fd'
            x_used = [y; (y - y_prev)/c.Ts];
            xh = x_used;
        case 'filtered'
            w_filt = (1 - a_filt)*w_filt + a_filt*(y - y_prev)/c.Ts;
            x_used = [y; w_filt];
            xh = x_used;
    end
    y_prev = y;

    % ---- control law ------------------------------------------------
    u = -c.K * x_used;
    if c.use_sat, u = min(max(u, c.umin), c.umax); end
    u_plant = u + c.dist_std*randn;              % input disturbance v[k]

    % ---- log ----------------------------------------------------------
    x_k(:,k) = x;  xhat(:,k) = xh;  u_log(k) = u;

    % ---- plant: hold u constant for one sample (ZOH) ------------------
    [ts, xs] = ode45(@(t, xx) pendulum_ode(xx, u_plant, c.p), ...
                     [tk(k), tk(k) + c.Ts], x, opts);
    t_cell{k} = ts(1:end-1);
    x_cell{k} = xs(1:end-1, :);
    x = xs(end, :)';
    u_prev = u;                                  % Unit Delay

    if abs(x(1)) > 4*pi, failed = true; break; end   % clearly fell over
end
if failed                                        % trim unused samples
    x_k = x_k(:, 1:k);  xhat = xhat(:, 1:k);  u_log = u_log(1:k);  tk = tk(1:k);
end

% ---- performance metrics --------------------------------------------
t_all = vertcat(t_cell{:});  x_all = vertcat(x_cell{:});
r.t = t_all';  r.x = x_all';  r.tk = tk;  r.x_k = x_k;  r.xhat = xhat;  r.u = u_log;
r.max_theta = max(abs(x_k(1,:)));
r.max_u     = max(abs(u_log));
r.effort    = sum(u_log.^2) * c.Ts;
band = 0.01;                                     % settled if |theta| < 0.01 rad
idx  = find(abs(x_k(1,:)) > band, 1, 'last');
if isempty(idx), r.t_settle = 0;
elseif idx == numel(tk), r.t_settle = Inf;
else, r.t_settle = tk(idx + 1);
end
r.success = ~failed && all(abs(x_k(1,:)) < pi/2) && abs(x_k(1,end)) < 0.05;
late = tk > 2;                                   % "steady state" part
r.rms_theta_end = sqrt(mean(x_k(1, late).^2));
r.rms_u_end     = sqrt(mean(u_log(late).^2));
r.rms_est_err   = sqrt(mean((x_k(2,:) - xhat(2,:)).^2));
end

function plot_run(r, ttl)
% Standard 3-panel plot: angle, velocity (true + estimate), control
figure('Name', ttl);
subplot(3,1,1); plot(r.t, r.x(1,:), 'LineWidth', 1.2); hold on;
plot(r.tk, r.xhat(1,:), '--'); grid on; ylabel('\theta [rad]');
legend('true', 'used by controller');
subplot(3,1,2); plot(r.t, r.x(2,:), 'LineWidth', 1.2); hold on;
plot(r.tk, r.xhat(2,:), '--'); grid on; ylabel('\omega [rad/s]');
legend('true', 'used by controller');
subplot(3,1,3); stairs(r.tk, r.u, 'LineWidth', 1.2); grid on;
ylabel('u [N m]'); xlabel('time [s]');
sgtitle(ttl);
end

function print_metrics(name, r)
fprintf('[%s] stable=%d  settle=%.2f s  max|theta|=%.3f rad  max|u|=%.2f N m  effort=%.3f\n', ...
        name, r.success, r.t_settle, r.max_theta, r.max_u, r.effort);
end

function save_fig(name)
file = fullfile('figures', [name '.png']);
try
    exportgraphics(gcf, file, 'Resolution', 150);   % MATLAB R2020a+
catch
    saveas(gcf, file);                              % older versions
end
end
