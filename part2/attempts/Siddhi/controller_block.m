function u  = controller(x,t,K_lqr,Fmax, Ts, Knom, Unom, Xnom, P_lqr, c_cap)
% Code of the MATLAB Function block in cartpole.slx (this file is a copy
% for reading / pasting; the model already contains it).
%
% x = [p; v; theta; omega]
% Angular measurements are in rad and rad/s. You need to implement a angle
% wrapping for theta to consider periodicity. You can use atan2 for this
% purpose.
% t is the time provided by the digital clock and Ts is the sampling time
% K_lqr is the infinite horizon LQR gain for the upright position
% Fmax is the maximum force
% Knom, Unom, Xnom are nominal gains, controls, and states determined by
% iLQR algorithm
% P_lqr, c_cap: capture criterion  x_w' * P_lqr * x_w <= c_cap
%   (P_lqr = Riccati matrix of the upright LQR, c_cap from the capture-zone
%    experiment, x_w = state with wrapped angle)
%
% ctrlMode: 0 = swing-up (track iLQR trajectory)
%           1 = upright LQR (latched once the capture criterion holds)
%           2 = failure (trajectory ended without capture) -> u = 0
%
% Angle wrapping:
%  - tracking: the error theta - theta_nom is wrapped to (-pi, pi], so a
%    measured angle that differs from the nominal one by 2*pi is treated
%    as the same configuration (the nominal angle itself is continuous,
%    pi -> 0).
%  - capture check and LQR: theta is wrapped to (-pi, pi], so upright is
%    always theta = 0, also when the plant reports 2*pi or -2*pi.

persistent ctrlMode
if isempty(ctrlMode)
    ctrlMode = 0;
end

N = numel(Unom);

% wrapped state for capture check / LQR
xw    = x(:);
xw(3) = atan2(sin(x(3)), cos(x(3)));

% sample index along the nominal trajectory (k = 1 <-> t = 0)
k = floor(t/Ts + 0.5) + 1;

if ctrlMode == 0
    if xw'*P_lqr*xw <= c_cap
        ctrlMode = 1;                                   % capture -> LQR
        fprintf('Capture at t = %.2f s -> switching to upright LQR\n', t);
    elseif k > N
        ctrlMode = 2;                                   % trajectory exhausted
        fprintf('SWING-UP FAILED: trajectory ended at t = %.2f s without reaching the capture zone (V = %.1f > %.1f)\n', ...
                t, xw'*P_lqr*xw, c_cap);
    end
end

if ctrlMode == 0
    k  = min(max(k,1),N);
    e    = x(:) - Xnom(:,k);
    e(3) = atan2(sin(e(3)), cos(e(3)));             % wrapped angle error
    u    = Unom(k) - Knom(:,:,k)*e;                 % u = ubar - K (x - xbar)
elseif ctrlMode == 1
    u = -K_lqr*xw;                                  % upright LQR
else
    u = 0;                                          % failed: no force
end

u = min(Fmax, max(-Fmax, u));                       % force limit
end
