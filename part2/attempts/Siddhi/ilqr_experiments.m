%% iLQR experiments (Sec. 4.1)
% Run setup_cartpole.m first (for model, Ts, x0, Q, R, Qf, c_cap, P_lqr).
% Each trial reports whether the nominal trajectory is ACCEPTED:
%   terminal state in the capture zone (x_N' P x_N <= c_cap) AND
%   max |p| <= xmax AND |u| <= Fmax (always true: forces are clipped).
% Note: P_lqr is recomputed for each Fmax (R_lqr = 1/Fmax^2).

c_cap = 55;   % capture threshold (also if an LQR-only run set it to Inf)
base = struct('Fmax',20,'xmax',0.5,'T',2.0,'Q',diag([10 0.1 0.1 0.1]), ...
              'R',1e-3,'Qf',1000*diag([1 1 10 1]),'guess','zero','boxed',true);
x0s = [0;0;pi;0];

fprintf('\n== backward pass: limit-aware vs. plain recursion (clipping only)\n');
for b = [true false]
    s = base; s.boxed = b; trial(s, sprintf('boxed = %d',b), x0s, model, Ts, c_cap);
end

fprintf('\n== maximum force\n');
for F = [5 8 10 12 15 20 30]
    s = base; s.Fmax = F; trial(s, sprintf('Fmax = %g',F), x0s, model, Ts, c_cap);
end
for F = [5 8 10]
    s = base; s.Fmax = F; s.T = 4; trial(s, sprintf('Fmax = %g, T = 4',F), x0s, model, Ts, c_cap);
end

fprintf('\n== cart limit\n');
for c = {[0.5 10 2],[0.3 10 2],[0.2 10 2],[0.2 100 2],[0.2 1000 2],[0.15 1000 3],[0.15 5000 3],[0.1 5000 4]}
    v = c{1}; s = base; s.xmax = v(1); s.Q(1,1) = v(2); s.T = v(3);
    trial(s, sprintf('xmax = %.2f, Q_p = %g, T = %g',v), x0s, model, Ts, c_cap);
end

fprintf('\n== weights\n');
W = {'baseline',            base.Q, base.R, base.Qf;
     'R = 0.01',            base.Q, 1e-2,   base.Qf;
     'R = 0.1',             base.Q, 1e-1,   base.Qf;
     'Qf / 10',             base.Q, base.R, base.Qf/10;
     'Qf / 100',            base.Q, base.R, base.Qf/100;
     'Qf * 10',             base.Q, base.R, base.Qf*10;
     'Q_theta = 1',         diag([10 0.1 1 0.1]),  base.R, base.Qf;
     'Q_theta = 10',        diag([10 0.1 10 0.1]), base.R, base.Qf;
     'Q_p = 0',             diag([0 0.1 0.1 0.1]), base.R, base.Qf;
     'Qf_omega small',      base.Q, base.R, 1000*diag([1 1 10 0.01])};
for i = 1:size(W,1)
    s = base; s.Q = W{i,2}; s.R = W{i,3}; s.Qf = W{i,4};
    trial(s, W{i,1}, x0s, model, Ts, c_cap);
end

fprintf('\n== horizon\n');
for T = [0.8 1.0 1.2 1.5 2 3 4 6]
    s = base; s.T = T; trial(s, sprintf('T = %g s',T), x0s, model, Ts, c_cap);
end

fprintf('\n== initial guess\n');
for gname = {'zero','sine','const','random'}
    for T = [2 3]
        s = base; s.guess = gname{1}; s.T = T;
        trial(s, sprintf('%s, T = %g',gname{1},T), x0s, model, Ts, c_cap);
    end
end

%% local functions
function r = trial(s, label, x0, model, Ts, c_cap)
N = round(s.T/Ts); t = (0:N-1)*Ts;
switch s.guess
    case 'zero',   U0 = zeros(1,N);
    case 'sine',   U0 = 0.5*s.Fmax*sin(2*pi*t/1.0);   % pumping at ~ natural freq.
    case 'const',  U0 = 0.5*s.Fmax*ones(1,N);
    case 'random', rng(1); U0 = s.Fmax*(2*rand(1,N)-1);
end
% LQR / capture zone for this force limit
[A,B] = upright_lin(model); Mx = expm([A B; zeros(1,5)]*Ts);
if exist('dlqr','file')
    [~,P] = dlqr(Mx(1:4,1:4), Mx(1:4,5), diag([4 1 25 0.25]), 1/s.Fmax^2);
else
    [~,P] = dlqr_iter(Mx(1:4,1:4), Mx(1:4,5), diag([4 1 25 0.25]), 1/s.Fmax^2);
end

opts = struct('maxIter',200,'tol',1e-4,'boxed',s.boxed);
tic; [X,U,~,Jh,info] = ilqr(x0,U0,model,Ts,s.Q,s.R,s.Qf,s.Fmax,opts); el = toc;
xN = X(:,end); xN(3) = atan2(sin(xN(3)),cos(xN(3)));
V = xN'*P*xN; pm = max(abs(X(1,:)));
why = '';
if V > c_cap,     why = [why ' terminal state outside capture zone;']; end
if pm > s.xmax,   why = [why ' cart limit exceeded;']; end
r.ok = isempty(why);
fprintf('%-30s %-8s %-18s it=%3d J=%9.1f V_N=%8.2f max|p|=%.3f thN=%+.2f omN=%+.2f (%.1fs)%s\n', ...
    label, ternary(r.ok,'ACCEPT','REJECT'), info.status, info.iterations, Jh(end), V, pm, ...
    X(3,end), X(4,end), el, why);
end

function [A,B] = upright_lin(p)
J = p.m*p.ell^2; D = (p.M+p.m)*J - (p.m*p.ell)^2;
A = [0 1 0 0; 0 -J*p.bc/D p.m^2*p.g*p.ell^2/D -p.m*p.ell*p.bp/D;
     0 0 0 1; 0 -p.m*p.ell*p.bc/D (p.M+p.m)*p.m*p.g*p.ell/D -(p.M+p.m)*p.bp/D];
B = [0; J/D; 0; p.m*p.ell/D];
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end

function [K,P] = dlqr_iter(A,B,Q,R)
P = Q;
for i = 1:100000
    K = (R + B'*P*B) \ (B'*P*A);
    Pn = Q + A'*P*(A - B*K);
    if norm(Pn - P, 'fro') < 1e-10*norm(P,'fro'), P = Pn; break, end
    P = Pn;
end
K = (R + B'*P*B) \ (B'*P*A);
end
