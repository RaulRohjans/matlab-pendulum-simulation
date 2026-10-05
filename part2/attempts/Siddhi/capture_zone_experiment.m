%% Empirical capture zone of the saturated upright LQR (Sec. 2.2)
% Simulates the sampled, saturated LQR u = -K_lqr x on the NONLINEAR model
% from many initial states and records which are stabilized.
% Success: after 6 s |theta|<0.02, |omega|<0.05, |p|<0.05, |v|<0.05
%          and |p| <= xmax during the whole transient (cart limit).
% Result -> capture criterion  x'P_lqr x <= c_cap  used in the controller.
%
% Run setup_cartpole.m first (or at least its LQR section) so that model,
% Ts, Fmax, xmax, K_lqr, P_lqr exist in the workspace.

c_cap = 55;   % capture threshold (also if an LQR-only run set it to Inf)
Tsim = 6;  nSteps = round(Tsim/Ts);

%% 1) 2-D slices (for plots)
th = linspace(-1.2,1.2,121);  om = linspace(-6,6,121);
[TH,OM] = meshgrid(th,om);
X0 = [zeros(2,numel(TH)); TH(:)'; OM(:)'];
okTO = reshape(sim_lqr(X0,K_lqr,model,Ts,Fmax,xmax,nSteps), size(TH));

pp = linspace(-0.5,0.5,101);  vv = linspace(-3,3,121);
[PP,VV] = meshgrid(pp,vv);
X0 = [PP(:)'; VV(:)'; zeros(2,numel(PP))];
okPV = reshape(sim_lqr(X0,K_lqr,model,Ts,Fmax,xmax,nSteps), size(PP));

%% 2) random 4-D sampling -> largest safe level set of x'P x
rng(0); n = 20000;
X0 = [ -0.4 + 0.8*rand(1,n);  -2 + 4*rand(1,n); ...
       -1.0 + 2.0*rand(1,n);  -5 + 10*rand(1,n)];
ok = sim_lqr(X0,K_lqr,model,Ts,Fmax,xmax,nSteps);
V  = sum(X0.*(P_lqr*X0),1);
c_star = min(V(~ok));
fprintf('Sampled %d states, %.1f%% stabilized.\n', n, 100*mean(ok));
fprintf('Largest level set with no failure: c* = %.2f\n', c_star);
for c = [0.5 0.75 1]*c_star
    in = V <= c;
    fprintf('  c = %6.2f : %5d samples inside, success rate %.4f\n', c, nnz(in), mean(ok(in)));
end
fprintf('Chosen c_cap = 0.75*c* = %.1f\n', 0.75*c_star);

% along the axes (others zero)
names = {'p','v','theta','omega'}; rngs = {0:0.01:0.5, 0:0.05:3, 0:0.01:1.5, 0:0.05:8};
for i = 1:4
    r = rngs{i}; Xa = zeros(4,numel(r)); Xa(i,:) = r;
    oka = sim_lqr(Xa,K_lqr,model,Ts,Fmax,xmax,nSteps);
    j = find(~oka,1); if isempty(j), lim = r(end); else, lim = r(max(j-1,1)); end
    fprintf('largest stabilized |%s| along its axis: %.2f\n', names{i}, lim);
end

%% 3) plots with the ellipse x'P x = c_cap
figure('Name','Capture zone theta-omega (p=v=0)');
imagesc(th,om,double(okTO)); axis xy; colormap([0.85 0.85 0.85; 0.4 0.7 0.4]); hold on
Pto = P_lqr([3 4],[3 4]);
fcontour(@(a,b) Pto(1,1)*a.^2 + 2*Pto(1,2)*a.*b + Pto(2,2)*b.^2, ...
         [th([1 end]) om([1 end])], 'LevelList', c_cap, 'LineColor','k','LineWidth',1.5);
xlabel('\theta [rad]'); ylabel('\omega [rad/s]');
title('green: stabilized by saturated LQR; black: x^TPx = c_{cap}');

figure('Name','Capture zone p-v (theta=omega=0)');
imagesc(pp,vv,double(okPV)); axis xy; colormap([0.85 0.85 0.85; 0.4 0.7 0.4]); hold on
Ppv = P_lqr([1 2],[1 2]);
fcontour(@(a,b) Ppv(1,1)*a.^2 + 2*Ppv(1,2)*a.*b + Ppv(2,2)*b.^2, ...
         [pp([1 end]) vv([1 end])], 'LevelList', c_cap, 'LineColor','k','LineWidth',1.5);
xlabel('p [m]'); ylabel('v [m/s]');
title('green: stabilized by saturated LQR; black: x^TPx = c_{cap}');

%% local functions
function ok = sim_lqr(X0,K,model,Ts,Fmax,xmax,nSteps)
% vectorized over columns of X0 (same equations as cartpole_dynamics + RK4)
x = X0; pmax = abs(x(1,:));
for k = 1:nSteps
    xw = x; xw(3,:) = atan2(sin(x(3,:)),cos(x(3,:)));
    u  = min(Fmax,max(-Fmax,-K*xw));
    a = f(x,u,model); b = f(x+Ts*a/2,u,model);
    c = f(x+Ts*b/2,u,model); d = f(x+Ts*c,u,model);
    x = x + Ts*(a+2*b+2*c+d)/6;
    pmax = max(pmax,abs(x(1,:)));
end
thw = atan2(sin(x(3,:)),cos(x(3,:)));
ok = abs(thw)<0.02 & abs(x(4,:))<0.05 & abs(x(1,:))<0.05 & abs(x(2,:))<0.05 & pmax<=xmax;
end

function dx = f(x,u,p)
v = x(2,:); th = x(3,:); om = x(4,:);
a11 = p.M+p.m; a12 = -p.m*p.ell*cos(th); a22 = p.m*p.ell^2;
r1 = u - p.bc*v - p.m*p.ell*sin(th).*om.^2;
r2 = p.m*p.g*p.ell*sin(th) - p.bp*om;
det = a11*a22 - a12.^2;
dx = [v; (a22*r1 - a12.*r2)./det; om; (-a12.*r1 + a11*r2)./det];
end
