function [success,summary,X] = simulate_saturated_lqr(X0,K,model,Ts,Fmax,xmax,Tsim)
%SIMULATE_SATURATED_LQR Vectorized nonlinear rollouts for capture studies.
% Each column of X0 is an initial state. A rollout succeeds when the cart
% respects xmax and the final state has settled close to upright.

if nargin < 8 || isempty(Tsim)
    Tsim = 6;
end
if size(X0,1) ~= 4
    error('CartPole:InitialStateShape','X0 must have four rows.');
end

X = X0;
maximumCartPosition = abs(X(1,:));
finiteTrajectory = all(isfinite(X),1);
nSteps = round(Tsim/Ts);

for k = 1:nSteps
    wrapped = X;
    wrapped(3,:) = atan2(sin(X(3,:)),cos(X(3,:)));
    force = min(Fmax,max(-Fmax,-K*wrapped));

    a = vector_dynamics(X,force,model);
    b = vector_dynamics(X+Ts*a/2,force,model);
    c = vector_dynamics(X+Ts*b/2,force,model);
    d = vector_dynamics(X+Ts*c,force,model);
    X = X + Ts*(a+2*b+2*c+d)/6;

    maximumCartPosition = max(maximumCartPosition,abs(X(1,:)));
    finiteTrajectory = finiteTrajectory & all(isfinite(X),1);
end

theta = atan2(sin(X(3,:)),cos(X(3,:)));
success = finiteTrajectory & maximumCartPosition <= xmax ...
    & abs(X(1,:)) < 0.05 & abs(X(2,:)) < 0.05 ...
    & abs(theta) < deg2rad(2) & abs(X(4,:)) < 0.08;

summary.maximumCartPosition = maximumCartPosition;
summary.finalWrappedAngle = theta;
summary.finalVelocity = X(2,:);
summary.finalAngularVelocity = X(4,:);
end

function dx = vector_dynamics(x,u,p)
v = x(2,:);
theta = x(3,:);
omega = x(4,:);
a11 = p.M+p.m;
a12 = -p.m*p.ell*cos(theta);
a22 = p.m*p.ell^2;
r1 = u-p.bc*v-p.m*p.ell*sin(theta).*omega.^2;
r2 = p.m*p.g*p.ell*sin(theta)-p.bp*omega;
determinant = a11*a22-a12.^2;
dx = [v;
      (a22*r1-a12.*r2)./determinant;
      omega;
      (-a12.*r1+a11*r2)./determinant];
end
