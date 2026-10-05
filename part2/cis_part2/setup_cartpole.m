%% Cart-pole swing-up (iLQR) + upright stabilization (LQR)
% Run this script. It designs the LQR, optimizes the swing-up with iLQR,
% checks the nominal trajectory and simulates the Simscape model.
%
% run_mode = 'swingup'  : full task (start hanging down, iLQR + LQR)
% run_mode = 'lqr_only' : Sec. 2.1 validation, start near upright and use
%                         only the LQR (capture criterion disabled)
clear; close all; clc
run_mode = 'swingup';

%% Physical parameters used by Simscape
%%%%% You can vary the parameters %%%%%%
param.m  = 1.0;     % Pendulum mass [kg]
param.M  = 1.0;     % Cart mass [kg]
param.L  = 0.5;     % Full pole length [m]
param.bc = 0.1;     % Cart damping [N*s/m]
param.bp = 0.05;    % Pivot damping [N*m*s/rad]
g=9.81;

%% Initial conditions: hanging down, stationary
%%%%% Vary the initial state %%%%%%
param.p0     = 0;
param.v0     = 0;
param.theta0 = 180; % Simscape joint target is in degrees (0 upright, 180 down)
param.omega0 = 0;
if strcmp(run_mode,'lqr_only')
    param.theta0 = 15;   % e.g. 15 deg from upright (inside the capture zone)
end

x0 = [param.p0; param.v0;
      deg2rad(param.theta0);param.omega0];

%% Parameters for nonlinear dynamics
model.M   = param.M;
model.m   = param.m;
model.ell = param.L/2;     % Pivot-to-point-mass distance [m]
model.bc  = param.bc;
model.bp  = param.bp;
model.g   = g;

%%%%% You can vary these parameters %%%%%%
Ts   = 0.01;               % Sample time for discrete-time dynamics [s]
Fmax = 20;                 % Force saturation [N]
xmax = 0.5;                % Cart-travel limit [m] (can be arbitrarily large)

%% Continuous-time linearization at upright
% Upright = 0; positive angles are counterclockwise from 0
%%% You do not need to change this part %%%

M   = model.M;
m   = model.m;
ell = model.ell;

J = m*ell^2;
Delta = (M+m)*J - (m*ell)^2;

A = [0, 1, 0, 0;
     0, -J*model.bc/Delta, ...
         m^2*g*ell^2/Delta, -m*ell*model.bp/Delta;
     0, 0, 0, 1;
     0, -m*ell*model.bc/Delta, ...
         (M+m)*m*g*ell/Delta, -(M+m)*model.bp/Delta];

B = [0;
     J/Delta;
     0;
     m*ell/Delta];

%% Discrete-time LQR
%%% Obtain a discrete-time model of the continuous time system
% Zero-order hold (input constant between samples):
%   Ad = e^(A Ts),  Bd = int_0^Ts e^(A s) ds B
% computed exactly with one matrix exponential of the augmented matrix
% [A B; 0 0]. Equivalent to c2d(ss(A,B,eye(4),0),Ts,'zoh').
Maug = expm([A, B; zeros(1,5)]*Ts);
Ad = Maug(1:4,1:4);
Bd = Maug(1:4,5);

%%% Compute the discrete LQR gain, select and justify weight matrices
% Bryson's rule: Q_ii = 1/(max acceptable x_i)^2, R = 1/(max acceptable u)^2
%   p     : 0.5 m     (= cart-travel limit xmax)
%   v     : 1.0 m/s
%   theta : 0.2 rad   (~11 deg, the quantity that must stay small)
%   omega : 2.0 rad/s
%   u     : Fmax      (use the whole force range before saturating)
Q_lqr = diag([1/0.5^2, 1/1^2, 1/0.2^2, 1/2^2]);
R_lqr = 1/Fmax^2;

% K_lqr: u = -K_lqr x ; P_lqr: Riccati solution (cost-to-go x'P x),
% used for the capture criterion  x'P_lqr x <= c_cap
if exist('dlqr','file')                       % Control System Toolbox
    [K_lqr, P_lqr] = dlqr(Ad, Bd, Q_lqr, R_lqr);
else                                          % fallback: Riccati iteration
    [K_lqr, P_lqr] = dlqr_iter(Ad, Bd, Q_lqr, R_lqr);
end
fprintf('K_lqr = [%s]\n', num2str(K_lqr,'%9.3f'));
fprintf('closed-loop |eig(Ad-Bd*K)| max = %.4f\n', max(abs(eig(Ad-Bd*K_lqr))));

% Capture criterion from capture_zone_experiment.m: with saturated LQR and
% |p| <= xmax, every sampled state with x'P x <= c_star = 73.1 was
% stabilized. A 25% margin is used for model mismatch -> c_cap = 55.
c_cap = 55;

%% iLQR setup
if strcmp(run_mode,'swingup')
T = 2.0;           % Swing-up horizon [s]
N = round(T/Ts);
T = N*Ts;

% Initial guess
U0 = zeros(1,N);   % no force; the linearization at the hanging state
                   % still gives a nonzero gradient. A pumping guess
                   % 0.5*Fmax*sin(2*pi*t) converges to a similar one
                   % (see ilqr_experiments.m).

% Swing-up weight matrices
%  - small running weights on v, theta, omega: the pendulum has to move far
%    from the origin, so they should not fight the swing-up
%  - Q_p = 10 keeps the cart near the centre (cart limit)
%  - large terminal weight, mostly on theta: end near upright and slow
%  - small R: the force limit, not R, is what bounds the force
Q  = diag([10, 0.1, 0.1, 0.1]);
R  = 1e-3;
Qf = 1000*diag([1, 1, 10, 1]);

% Force margin: the nominal swing-up is planned with 90% of the force
% limit, so 2 N remain for the feedback -K_k (x - xbar) when the real
% plant deviates from the model. Without the margin the nominal force is
% at +-Fmax most of the time and the feedback saturates (see report).
Fplan = 0.9*Fmax;

%% Optimize swing-up
%%% Implement ilqr function %%%
opts = struct('maxIter',200,'tol',1e-4,'verbose',true);
tic
[Xnom,Unom,Knom,costHistory,info] = ilqr(x0,U0,model,Ts,Q,R,Qf,Fplan,opts);
fprintf('iLQR: %s, %d iterations, J = %.2f, %.1f s\n', ...
        info.status, info.iterations, costHistory(end), toc);

%% Check nominal trajectory

%%% Check if the nominal trajectory Xnum, Unom returned by iLQR satisfies
%%% that the terminal state is in the capture zone, the cart displacement
%%% is within limits, the force is bounded (it should be)
xN    = Xnom(:,end);
xN(3) = atan2(sin(xN(3)),cos(xN(3)));
V_N   = xN'*P_lqr*xN;
checks.capture = V_N <= c_cap;
checks.cart    = max(abs(Xnom(1,:))) <= xmax;
checks.force   = max(abs(Unom)) <= Fmax + 1e-9;
fprintf('Terminal state: p=%.3f v=%.3f theta=%.3f omega=%.3f, V_N = %.2f (c_cap = %g)\n', ...
        xN, V_N, c_cap);
fprintf('max|p| = %.3f m (xmax = %g), max|u| = %.2f N (Fplan = %g, Fmax = %g)\n', ...
        max(abs(Xnom(1,:))), xmax, max(abs(Unom)), Fplan, Fmax);
if ~checks.capture
    warning(['Nominal trajectory REJECTED: terminal state outside capture zone. ' ...
             'Try a longer horizon T, a larger Qf (theta/omega entries) or more force.']);
end
if ~checks.cart
    warning(['Nominal trajectory REJECTED: cart limit exceeded. ' ...
             'Increase Q(1,1) (cart position weight) and/or the horizon T.']);
end
if ~checks.force
    warning('Nominal trajectory REJECTED: force limit violated.');
end
if checks.capture && checks.cart && checks.force
    fprintf('Nominal trajectory ACCEPTED.\n');
end

%% Plot results

%%%% You can plot the nominal state and control trajectories. These can be
%%%% informative for changing iLQR parameters for finding a feasible
%%%% trajectory.
tn = (0:N)*Ts;
figure('Name','iLQR nominal trajectory');
lbl = {'p [m]','v [m/s]','\theta [rad]','\omega [rad/s]'};
for i = 1:4
    subplot(5,1,i); plot(tn, Xnom(i,:), 'LineWidth',1.2); grid on; ylabel(lbl{i});
end
subplot(5,1,1); hold on; yline([-xmax xmax],'r--');
subplot(5,1,5); stairs(tn(1:end-1), Unom, 'LineWidth',1.2); grid on; hold on;
yline([-Fmax Fmax],'r--'); yline([-Fplan Fplan],'k:'); ylabel('u [N]'); xlabel('t [s]');

figure('Name','iLQR cost');
semilogy(0:numel(costHistory)-1, costHistory, 'o-'); grid on;
xlabel('iteration'); ylabel('J');

else
    % LQR-only test: dummy nominal trajectory, always use the LQR
    Xnom = zeros(4,2); Unom = zeros(1,2); Knom = zeros(1,4,2);   % 3-D like the real gains
    c_cap = Inf;
end

%% Simulation
out=sim("cartpole");

%%% Plot additional results. For example, comparing nominal and resulting
%%% swing-up trajectory.
ts = out.x.Time;
xs = squeeze(out.x.Data);
if size(xs,1) == 4 && size(xs,2) ~= 4, xs = xs.'; end   % -> (samples x 4)
us = squeeze(out.u.Data);
tu = out.u.Time;

% reconstruct when the controller switched (same rule as in the block)
th_w = atan2(sin(xs(:,3)),cos(xs(:,3)));
xw   = [xs(:,1:2), th_w, xs(:,4)];
Vs   = sum((xw*P_lqr).*xw, 2);
iCap = find(Vs <= c_cap, 1);
if isempty(iCap)
    fprintf(2,'FAILURE: the state never entered the capture zone.\n');
    tCap = NaN;
else
    tCap = ts(iCap);
    fprintf('Captured at t = %.2f s.\n', tCap);
end
ok = abs(th_w(end)) < 0.05 && abs(xs(end,4)) < 0.1;
fprintf('Final state: p=%.3f v=%.3f theta(wrapped)=%.4f omega=%.4f -> %s\n', ...
        xs(end,1), xs(end,2), th_w(end), xs(end,4), ...
        ternary(ok,'STABILIZED upright','NOT stabilized'));
fprintf('max|p| = %.3f m (xmax = %g)\n', max(abs(xs(:,1))), xmax);

figure('Name','Simscape closed loop');
for i = 1:4
    subplot(5,1,i);
    if i == 3, plot(ts, unwrap(xs(:,i)), 'LineWidth',1.2);
    else,      plot(ts, xs(:,i), 'LineWidth',1.2); end
    hold on; grid on; ylabel(lbl_or(i));
    if strcmp(run_mode,'swingup'), plot((0:numel(Unom))*Ts, Xnom(i,:), '--'); end
    if ~isnan(tCap), xline(tCap,'k:'); end
end
subplot(5,1,1); legend('Simscape','iLQR nominal','capture','Location','best');
subplot(5,1,5); stairs(tu, us, 'LineWidth',1.2); grid on; hold on;
if strcmp(run_mode,'swingup'), stairs((0:numel(Unom)-1)*Ts, Unom, '--'); end
ylabel('u [N]'); xlabel('t [s]');
if ~isnan(tCap), xline(tCap,'k:'); end

% zoom on the swing-up phase
if strcmp(run_mode,'swingup'), for i = 1:5, subplot(5,1,i); xlim([0 max(4,T+1)]); end, end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
function s = lbl_or(i)
L = {'p [m]','v [m/s]','\theta [rad]','\omega [rad/s]'}; s = L{i};
end
function [K,P] = dlqr_iter(A,B,Q,R)
% Discrete algebraic Riccati equation by fixed-point iteration
P = Q;
for i = 1:100000
    K = (R + B'*P*B) \ (B'*P*A);
    Pn = Q + A'*P*(A - B*K);
    if norm(Pn - P, 'fro') < 1e-10*norm(P,'fro'), P = Pn; break, end
    P = Pn;
end
K = (R + B'*P*B) \ (B'*P*A);
end
